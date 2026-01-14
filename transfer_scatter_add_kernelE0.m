function [S_NL_E, diag] = transfer_scatter_add_kernelE0(kVals, edges, E, weight, centroidX, centroidY)

%function [S_NL_E, diag] = FV_energy_transfer_uniquePQ_exact(kVals, edges, E, weight, centroidX, centroidY)

kVals = double(kVals(:));
edges = double(edges(:));
E     = double(E(:));
weight    = double(weight);
centroidX = double(centroidX);
centroidY = double(centroidY);

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

diag.numTriadsUsed = 0;
diag.maxTriadEnergyResidual = 0;

for kj = 1:N
    kstar = krep(kj);

    for pj = 1:N
        for qj = pj:N

            wk = weight(pj,qj,kj);
            if wk <= 0, continue; end

            dv = wk * dk(kj);

            % ---- (pj,qj) contribution
            [dEk,dEp,dEq] = triad_energy_increment(kstar, pj, qj, kj, dv, krep, edges, centroidX, centroidY, E0_at);

            % Kahan add
            [dE(kj), cE(kj)] = kahan_add(dE(kj), cE(kj), dEk);
            [dE(pj), cE(pj)] = kahan_add(dE(pj), cE(pj), dEp);
            [dE(qj), cE(qj)] = kahan_add(dE(qj), cE(qj), dEq);

            diag.maxTriadEnergyResidual = max(diag.maxTriadEnergyResidual, abs(dEk+dEp+dEq));
            diag.numTriadsUsed = diag.numTriadsUsed + 1;

            % ---- add the swapped permutation explicitly if qj>pj
            if qj > pj
                % swapped: (qj,pj) uses centroidX(qj,pj,kj)=qstar, centroidY(...)=pstar if you mirrored correctly
                wk2 = weight(qj,pj,kj);
                if wk2 > 0
                    dv2 = wk2 * dk(kj);

                    [dEk2,dEq2,dEp2] = triad_energy_increment(kstar, qj, pj, kj, dv2, krep, edges, centroidX, centroidY, E0_at);
                    % Note: because we swapped p<->q indices, what comes back as "p-leg"
                    % must be added to bin qj, and "q-leg" to bin pj, hence the naming.

                    [dE(kj), cE(kj)] = kahan_add(dE(kj), cE(kj), dEk2);
                    [dE(pj), cE(pj)] = kahan_add(dE(pj), cE(pj), dEp2);
                    [dE(qj), cE(qj)] = kahan_add(dE(qj), cE(qj), dEq2);

                    diag.maxTriadEnergyResidual = max(diag.maxTriadEnergyResidual, abs(dEk2+dEq2+dEp2));
                    diag.numTriadsUsed = diag.numTriadsUsed + 1;
                end
            end
        end
    end
end

S_NL_E = dE ./ dk;
diag.FV_total_energy_transfer = sum(S_NL_E .* dk);   % = sum(dE)

end

% ------------------------------------------------------------
% Triad energy increment for one ordered (pj,qj,kj)
% Returns ENERGY increments (already multiplied by dv and Jacobians).
% ------------------------------------------------------------
function [dEk,dEp,dEq] = triad_energy_increment(kstar, pj, qj, kj, dv, krep, edges, centroidX, centroidY, E0_at)

% One consistent triplet point
pstar = centroidX(pj,qj,kj);
qstar = centroidY(pj,qj,kj);

if ~(isfinite(pstar) && isreal(pstar) && pstar>0)
    pstar = sqrt(edges(pj)*edges(pj+1));
end
if ~(isfinite(qstar) && isreal(qstar) && qstar>0)
    qstar = sqrt(edges(qj)*edges(qj+1));
end

E0k = E0_at(kstar);
E0p = E0_at(pstar);
E0q = E0_at(qstar);

% Your raw 3-leg kernel in E0
Sk_raw = 0.5*( E0q*(E0p-E0k) + E0p*(E0q-E0k) );
Sq_raw = 0.5*( E0k*(E0p-E0q) + E0p*(E0k-E0q) );
Sp_raw = 0.5*( E0k*(E0q-E0p) + E0q*(E0k-E0p) );

% Jacobians at same evaluation points
Jk = 4*pi*kstar^2;
Jp = 4*pi*pstar^2;
Jq = 4*pi*qstar^2;

% Energy-conserving delta
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
