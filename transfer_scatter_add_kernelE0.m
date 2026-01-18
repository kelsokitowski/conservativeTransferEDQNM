function [S_NL_E, FV_total_energy_transfer, maxTriadEnergyResidual, numTriadsUsed] = transfer_scatter_add_kernelE0(kVals, edges, E, weight_k, centroidX_k, centroidY_k, centroidZ_k, weight_p, centroidX_p, centroidY_p, weight_q, centroidX_q, centroidY_q, mu1, nu, t)
% transfer_scatter_add_kernelE0 (K-SLICE INTEGRATION)
%
% Integrates over k-slices only: each triad integrated once with proper volume weighting.
% The cyclic symmetry of dv (verified externally) ensures geometric correctness.
% Delta correction at each triad enforces local energy conservation.
%
% EDQNM Kernel: theta * 16*pi^2 * p^2*k^2*q * (xy+z^3) * E0q*(E0p-E0k)
% where theta includes relaxation time dynamics
%
% PARALLELIZATION STRATEGY (OpenMP/MPI):
%
% Each kj can be computed independently. Two approaches:
%
% APPROACH 1 - Thread-local arrays (best for OpenMP):
%   Allocate dE_local(kLength) per thread, accumulate locally with Kahan,
%   then combine all thread-local arrays at the end with Kahan summation.
%   This preserves full precision.
%
%   !$OMP PARALLEL PRIVATE(pj,qj,wk,dv_local,pstar,qstar,dEk,dEp,dEq, &
%   !$OMP                  dE_local,cE_local)
%   allocate(dE_local(kLength), cE_local(kLength))
%   dE_local = 0.0; cE_local = 0.0
%   !$OMP DO SCHEDULE(dynamic)
%   do kj = 1, kLength
%       ... compute contributions, add to dE_local using kahan_add ...
%   end do
%   !$OMP END DO
%   !$OMP CRITICAL
%   do kj = 1, kLength
%       call kahan_add(dE(kj), cE(kj), dE_local(kj), ...)
%   end do
%   !$OMP END CRITICAL
%   deallocate(dE_local, cE_local)
%   !$OMP END PARALLEL
%
% APPROACH 2 - MPI domain decomposition (best for distributed memory):
%   Each MPI rank computes S_NL_E for a range of kj indices:
%     rank r handles kj = kj_start(r) to kj_end(r)
%   Each rank returns its local S_NL_E(kj_start:kj_end).
%   No global reduction needed - each rank owns its kj bins.
%   Near-linear scaling up to kLength ranks.
%
% Outputs:
%   S_NL_E                    - Energy transfer rate (kLength x 1)
%   FV_total_energy_transfer  - Sum of all energy transfers (scalar)
%   maxTriadEnergyResidual    - Maximum residual across all triads (scalar)
%   numTriadsUsed             - Number of triads processed (scalar)

kVals = double(kVals(:));
edges = double(edges(:));
E     = double(E(:));
mu1   = double(mu1(:));
nu    = double(nu);
t     = double(t);

weight_k    = double(weight_k);
centroidX_k = double(centroidX_k);
centroidY_k = double(centroidY_k);

weight_p    = double(weight_p);
centroidX_p = double(centroidX_p);
centroidY_p = double(centroidY_p);

weight_q    = double(weight_q);
centroidX_q = double(centroidX_q);
centroidY_q = double(centroidY_q);

kLength = length(kVals);
dk = diff(edges(:));
krep = sqrt(edges(1:end-1).*edges(2:end));  % FV reps

% Precompute E0 and mu1 interpolation data (read-only, thread-safe)
E0  = E ./ (4*pi*kVals.^2);
E0s = max(E0, realmin);
logk  = log(kVals);
logE0 = log(E0s);

% mu1 interpolation data (linear in log-log space)
mu1s = max(mu1, realmin);
logmu1 = log(mu1s);

% Energy increment accumulators
% Note: Kahan summation needed for precision when accumulating ~60k triads
dE = zeros(kLength,1);
cE = zeros(kLength,1);  % Kahan compensation

% Diagnostic scalars (thread-safe with OpenMP REDUCTION)
numTriadsUsed = 0;
maxTriadEnergyResidual = 0.0;

% ===========================================================================
% K-SLICE LOOP: Parallelizable over kj (see OpenMP directives in header)
% ===========================================================================
for kj = 1:kLength
    for pj = 1:kLength
        for qj = pj:kLength
            wk = weight_k(pj,qj,kj);
            if wk <= 0, continue; end

            dv_local = wk * dk(kj);

            % Evaluate kernel at 3D centroid (kstar, pstar, qstar) for exact Sk+Sp+Sq=0
            kstar = centroidZ_k(pj,qj,kj);  % k-centroid
            pstar = centroidX_k(pj,qj,kj);  % p-centroid
            qstar = centroidY_k(pj,qj,kj);  % q-centroid

            [dEk,dEp,dEq] = triad_energy_increment_direct(kstar, pstar, qstar, dv_local, kj, pj, qj, logk, logE0, logmu1, nu, t);

            % Scatter-add to bins using Kahan summation for precision
            [dE(kj), cE(kj)] = kahan_add(dE(kj), cE(kj), dEk);
            [dE(pj), cE(pj)] = kahan_add(dE(pj), cE(pj), dEp);
            [dE(qj), cE(qj)] = kahan_add(dE(qj), cE(qj), dEq);

            maxTriadEnergyResidual = max(maxTriadEnergyResidual, abs(dEk+dEp+dEq));
            numTriadsUsed = numTriadsUsed + 1;

            % Mirror p<->q if needed
            %if qj > pj
            %    wk2 = weight_k(qj,pj,kj);
            %    if wk2 > 0
            %        dv2 = wk2 * dk(kj);
            %        pstar2 = centroidX_k(qj,pj,kj);  % swapped
            %        qstar2 = centroidY_k(qj,pj,kj);  % swapped

            %        [dEk2,dEp2,dEq2] = triad_energy_increment_direct(kstar, pstar2, qstar2, dv2, edges, E0_at);

            %        [dE(kj), cE(kj)] = kahan_add(dE(kj), cE(kj), dEk2);
            %        [dE(pj), cE(pj)] = kahan_add(dE(pj), cE(pj), dEp2);
            %        [dE(qj), cE(qj)] = kahan_add(dE(qj), cE(qj), dEq2);

            %        diag.maxTriadEnergyResidual = max(diag.maxTriadEnergyResidual, abs(dEk2+dEp2+dEq2));
            %        diag.numTriadsUsed = diag.numTriadsUsed + 1;
            %    end
            %end
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
FV_total_energy_transfer = sum(dE);  % Total energy change (should be ~0)

end

% ============================================================================
% HELPER FUNCTION: Triad energy increment
% ============================================================================
% Direct evaluation at (kstar, pstar, qstar)
% Returns ENERGY increments (already multiplied by dv and Jacobians).
function [dEk,dEp,dEq] = triad_energy_increment_direct(kstar, pstar, qstar, dv, kj, pj, qj, logk, logE0, logmu1, nu, t)

% Validate centroids
if ~(isfinite(pstar) && isreal(pstar) && pstar>0)
    error('Invalid pstar: %g', pstar);
end
if ~(isfinite(qstar) && isreal(qstar) && qstar>0)
    error('Invalid qstar: %g', qstar);
end
if ~(isfinite(kstar) && isreal(kstar) && kstar>0)
    error('Invalid kstar: %g', kstar);
end

% Interpolate E0 at triad evaluation points
E0k = interpolate_E0(kstar, logk, logE0);
E0p = interpolate_E0(pstar, logk, logE0);
E0q = interpolate_E0(qstar, logk, logE0);

% Interpolate mu1 at triad evaluation points
mu1_k = interpolate_mu1(kstar, logk, logmu1);
mu1_p = interpolate_mu1(pstar, logk, logmu1);
mu1_q = interpolate_mu1(qstar, logk, logmu1);

% Compute EDQNM kernel for each leg (two terms per leg)
% ALL arguments permuted: wavenumbers, spectral densities, mu1 values, indices
%
% k-leg: Sk_raw = 0.5*(kernel1(k,p,q) + kernel1(k,q,p))
% Cyclic (k,p,q): spectral from (p,q,k), wavenumbers (k,p,q), mu (k,p,q), indices (k,p,q)
term1_k = kernel1(E0p, E0q, E0k, kstar, pstar, qstar, mu1_k, mu1_p, mu1_q, kj, pj, qj, nu, t);
% Anticyclic (k,q,p): spectral from (q,p,k), wavenumbers (k,q,p), mu (k,q,p), indices (k,q,p)
term2_k = kernel1(E0q, E0p, E0k, kstar, qstar, pstar, mu1_k, mu1_q, mu1_p, kj, qj, pj, nu, t);
Sk_raw = 0.5 * (term1_k + term2_k);

% Diagnostic: print first few triads to check magnitudes
persistent triad_count;
if isempty(triad_count), triad_count = 0; end
triad_count = triad_count + 1;
if triad_count <= 5
    thetaVal_test = theta(nu, kstar, pstar, qstar, kj, pj, qj, t, mu1_k, mu1_p, mu1_q);
    damping_rate = nu*(kstar^2 + pstar^2 + qstar^2) + mu1_k + mu1_p + mu1_q;
    fprintf('Triad %d: kj=%d pj=%d qj=%d | k=%.3e p=%.3e q=%.3e\n', triad_count, kj, pj, qj, kstar, pstar, qstar);
    fprintf('  E0: k=%.3e p=%.3e q=%.3e | mu1: k=%.3e p=%.3e q=%.3e\n', E0k, E0p, E0q, mu1_k, mu1_p, mu1_q);
    fprintf('  theta=%.3e damping=%.3e t=%.3e nu=%.3e\n', thetaVal_test, damping_rate, t, nu);
    fprintf('  term1_k=%.3e term2_k=%.3e Sk_raw=%.3e | dv=%.3e\n', term1_k, term2_k, Sk_raw, dv);
end

% p-leg: Sp_raw = 0.5*(kernel1(p,q,k) + kernel1(p,k,q))
% Cyclic (p,q,k): spectral from (q,k,p), wavenumbers (p,q,k), mu (p,q,k), indices (p,q,k)
term1_p = kernel1(E0q, E0k, E0p, pstar, qstar, kstar, mu1_p, mu1_q, mu1_k, pj, qj, kj, nu, t);
% Anticyclic (p,k,q): spectral from (k,q,p), wavenumbers (p,k,q), mu (p,k,q), indices (p,k,q)
term2_p = kernel1(E0k, E0q, E0p, pstar, kstar, qstar, mu1_p, mu1_k, mu1_q, pj, kj, qj, nu, t);
Sp_raw = 0.5 * (term1_p + term2_p);

% q-leg: Sq_raw = 0.5*(kernel1(q,k,p) + kernel1(q,p,k))
% Cyclic (q,k,p): spectral from (k,p,q), wavenumbers (q,k,p), mu (q,k,p), indices (q,k,p)
term1_q = kernel1(E0k, E0p, E0q, qstar, kstar, pstar, mu1_q, mu1_k, mu1_p, qj, kj, pj, nu, t);
% Anticyclic (q,p,k): spectral from (p,k,q), wavenumbers (q,p,k), mu (q,p,k), indices (q,p,k)
term2_q = kernel1(E0p, E0k, E0q, qstar, pstar, kstar, mu1_q, mu1_p, mu1_k, qj, pj, kj, nu, t);
Sq_raw = 0.5 * (term1_q + term2_q);

% Jacobians at evaluation points
% Jk = 4*pi*kstar^2;
% Jp = 4*pi*pstar^2;
% Jq = 4*pi*qstar^2;
%Jk = 1; Jp = 1; Jq = 1;

% Energy-conserving delta correction
%delta = (Jk*Sk_raw + Jp*Sp_raw + Jq*Sq_raw) / (Jk + Jp + Jq);
delta = (Sk_raw + Sp_raw + Sq_raw) / (3.0);

Sk = Sk_raw - delta;
Sp = Sp_raw - delta;
Sq = Sq_raw - delta;

% Energy increments


dEk = Sk * dv;
dEp = Sp * dv;
dEq = Sq * dv;

% Final exact-zero projection (kills roundoff at triad level)
s = dEk + dEp + dEq;
dEk = dEk - s/3.0;
dEp = dEp - s/3.0;
dEq = dEq - s/3.0;

end

% ============================================================================
% HELPER FUNCTION: E0 interpolation (Fortran-compatible)
% ============================================================================
% Log-log linear interpolation of E0
% In Fortran: implement as a function that searches logk array and interpolates
function E0_val = interpolate_E0(k_eval, logk, logE0)
logk_eval = log(k_eval);
logE0_interp = interp1(logk, logE0, logk_eval, 'linear', 'extrap');
E0_val = exp(logE0_interp);
end

% ============================================================================
% HELPER FUNCTION: mu1 interpolation (Fortran-compatible)
% ============================================================================
% Log-log linear interpolation of mu1
% In Fortran: implement as a function that searches logk array and interpolates
function mu1_val = interpolate_mu1(k_eval, logk, logmu1)
logk_eval = log(k_eval);
logmu1_interp = interp1(logk, logmu1, logk_eval, 'linear', 'extrap');
mu1_val = exp(logmu1_interp);
end

% ============================================================================
% HELPER FUNCTION: EDQNM kernel (Fortran-compatible)
% ============================================================================
% Computes one term of EDQNM kernel: theta * 16*pi^2 * p^2*k^2*q * (xy+z^3) * E0_q*(E0_p-E0_k)
% where:
%   E0_p, E0_q, E0_k: spectral energy density at p, q, k (cyclic ordering convention)
%   k, p, q: triad wavenumbers (for geometry computation)
%   mu1_k, mu1_p, mu1_q: eddy damping at k, p, q
%   kj, pj, qj: bin indices
%   nu, t: viscosity and integration time
function kernel1Val = kernel1(E0_p, E0_q, E0_k, k, p, q, mu1_k, mu1_p, mu1_q, kj, pj, qj, nu, t)

% Triad geometry: cosines of angles (computed from k, p, q wavenumbers)
% x = cos(angle between k and p) = (k^2 + p^2 - q^2)/(2*k*p)
% y = cos(angle between k and q) = (k^2 + q^2 - p^2)/(2*k*q)
% z = cos(angle between p and q) = (p^2 + q^2 - k^2)/(2*p*q)
x = (k^2 + p^2 - q^2) / (2*k*p);
y = (k^2 + q^2 - p^2) / (2*k*q);
z = (p^2 + q^2 - k^2) / (2*p*q);

% Theta relaxation function (uses k, p, q wavenumbers)
thetaVal = theta(nu, k, p, q, kj, pj, qj, t, mu1_k, mu1_p, mu1_q);

% Full EDQNM kernel term
% Geometry: 16*pi^2 * p^2*k^2*q * (xy+z^3)
% Spectral: E0_q*(E0_p-E0_k)
kernel1Val = thetaVal * 16.0 * pi^2 * p^2 * k^2 * q * (x*y + z^3) * E0_q * (E0_p - E0_k);

end

% ============================================================================
% HELPER FUNCTION: Kahan summation
% ============================================================================
% Compensated summation for improved numerical accuracy
% See: Kahan, W. (1965). "Further remarks on reducing truncation errors"
function [sum_new, c_new] = kahan_add(sum_old, c_old, x)
y = x - c_old;
t = sum_old + y;
c_new = (t - sum_old) - y;
sum_new = t;
end

