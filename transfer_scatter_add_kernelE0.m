function [S_NL_E, FV_total_energy_transfer, maxTriadEnergyResidual, numTriadsUsed] = transfer_scatter_add_kernelE0(kVals, edges, E, weight_k, centroidX_k, centroidY_k, weight_p, centroidX_p, centroidY_p, weight_q, centroidX_q, centroidY_q)
% transfer_scatter_add_kernelE0 (K-SLICE INTEGRATION)
%
% Integrates over k-slices only: each triad integrated once with proper volume weighting.
% The cyclic symmetry of dv (verified externally) ensures geometric correctness.
% Delta correction at each triad enforces local energy conservation.
%
% Outputs:
%   S_NL_E                    - Energy transfer rate (N x 1)
%   FV_total_energy_transfer  - Sum of all energy transfers (scalar)
%   maxTriadEnergyResidual    - Maximum residual across all triads (scalar)
%   numTriadsUsed             - Number of triads processed (scalar)

kVals = double(kVals(:));
edges = double(edges(:));
E     = double(E(:));

weight_k    = double(weight_k);
centroidX_k = double(centroidX_k);
centroidY_k = double(centroidY_k);

weight_p    = double(weight_p);
centroidX_p = double(centroidX_p);
centroidY_p = double(centroidY_p);

weight_q    = double(weight_q);
centroidX_q = double(centroidX_q);
centroidY_q = double(centroidY_q);

N  = length(kVals);
dk = diff(edges(:));
krep = sqrt(edges(1:end-1).*edges(2:end));  % FV reps

% Kernel variable
E0  = E ./ (4*pi*kVals.^2);
E0s = max(E0, realmin);
logk  = log(kVals);
logE0 = log(E0s);
E0_at = @(x) exp(interp1(logk, logE0, log(x), 'linear', 'extrap'));

% Kahan accumulators for bin energy increments
dE = zeros(N,1);  cE = zeros(N,1);

% Diagnostic scalars
numTriadsUsed = 0;
maxTriadEnergyResidual = 0;

% ===========================================================================
% CONTRIBUTION 1: k-slices (integrate over p,q for fixed k)
% ===========================================================================
for kj = 1:N
    kstar = krep(kj);

    for pj = 1:N
        for qj = pj:N
            wk = weight_k(pj,qj,kj);
            if wk <= 0, continue; end

            dv = wk * dk(kj);

            % Evaluate kernel at (kstar, pstar, qstar) from k-slice centroids
            pstar = centroidX_k(pj,qj,kj);
            qstar = centroidY_k(pj,qj,kj);

            [dEk,dEp,dEq] = triad_energy_increment_direct(kstar, pstar, qstar, dv, edges, E0_at);

            % Scatter-add to bins (kj, pj, qj)
            [dE(kj), cE(kj)] = kahan_add(dE(kj), cE(kj), dEk);
            [dE(pj), cE(pj)] = kahan_add(dE(pj), cE(pj), dEp);
            [dE(qj), cE(qj)] = kahan_add(dE(qj), cE(qj), dEq);

            maxTriadEnergyResidual = max(maxTriadEnergyResidual, abs(dEk+dEp+dEq));
            numTriadsUsed = numTriadsUsed + 1;

            % Mirror p<->q if needed
            if qj > pj
                wk2 = weight_k(qj,pj,kj);
                if wk2 > 0
                    dv2 = wk2 * dk(kj);
                    pstar2 = centroidX_k(qj,pj,kj);  % swapped
                    qstar2 = centroidY_k(qj,pj,kj);  % swapped

                    [dEk2,dEp2,dEq2] = triad_energy_increment_direct(kstar, pstar2, qstar2, dv2, edges, E0_at);

                    [dE(kj), cE(kj)] = kahan_add(dE(kj), cE(kj), dEk2);
                    [dE(pj), cE(pj)] = kahan_add(dE(pj), cE(pj), dEp2);
                    [dE(qj), cE(qj)] = kahan_add(dE(qj), cE(qj), dEq2);

                    maxTriadEnergyResidual = max(maxTriadEnergyResidual, abs(dEk2+dEp2+dEq2));
                    numTriadsUsed = numTriadsUsed + 1;
                end
            end
        end
    end
end

% ===========================================================================
% NOTE: p-slices and q-slices commented out - using only k-slice integration
% ===========================================================================
% The cyclic symmetry of dv is verified for correctness, but we integrate
% each triad only once (via k-slice) to avoid triple-counting.
% This preserves physics: each triad integrated once with proper volume weighting.
%
% % CONTRIBUTION 2: p-slices (integrate over q,k for fixed p)
% % CONTRIBUTION 3: q-slices (integrate over k,p for fixed q)

S_NL_E = dE ./ dk;
FV_total_energy_transfer = sum(S_NL_E .* dk);   % = sum(dE)

end

% ------------------------------------------------------------
% Triad energy increment (direct evaluation at (kstar, pstar, qstar))
% Returns ENERGY increments (already multiplied by dv and Jacobians).
% ------------------------------------------------------------
function [dEk,dEp,dEq] = triad_energy_increment_direct(kstar, pstar, qstar, dv, edges, E0_at)

% Fallback to bin centers if centroids are invalid
if ~(isfinite(pstar) && isreal(pstar) && pstar>0)
    error('Invalid pstar: %g', pstar);
end
if ~(isfinite(qstar) && isreal(qstar) && qstar>0)
    error('Invalid qstar: %g', qstar);
end
if ~(isfinite(kstar) && isreal(kstar) && kstar>0)
    error('Invalid kstar: %g', kstar);
end

E0k = E0_at(kstar);
E0p = E0_at(pstar);
E0q = E0_at(qstar);

% Raw 3-leg kernel in E0 (toy kernel)
Sk_raw = 0.5*( E0q*(E0p-E0k) + E0p*(E0q-E0k) );
Sq_raw = 0.5*( E0k*(E0p-E0q) + E0p*(E0k-E0q) );
Sp_raw = 0.5*( E0k*(E0q-E0p) + E0q*(E0k-E0p) );

% Jacobians at evaluation points
Jk = 4*pi*kstar^2;
Jp = 4*pi*pstar^2;
Jq = 4*pi*qstar^2;

% Energy-conserving delta correction
delta = (Jk*Sk_raw + Jp*Sp_raw + Jq*Sq_raw) / (Jk + Jp + Jq);

Sk = Sk_raw - delta;
Sp = Sp_raw - delta;
Sq = Sq_raw - delta;

% Energy increments
dEk = Jk * Sk * dv;
dEp = Jp * Sp * dv;
dEq = Jq * Sq * dv;

% Final exact-zero projection (kills roundoff at triad level)
s = dEk + dEp + dEq;
dEk = dEk - s/3;
dEp = dEp - s/3;
dEq = dEq - s/3;

end

function [sum_new, c_new] = kahan_add(sum_old, c_old, x)
y = x - c_old;
t = sum_old + y;
c_new = (t - sum_old) - y;
sum_new = t;
end
