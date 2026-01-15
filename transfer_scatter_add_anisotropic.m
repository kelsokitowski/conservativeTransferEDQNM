function [S_NL_ISO, S_NL_DIR, S_NL_POL, ST_NL_ISO, ST_NL_DIR, SF_NL, diagnostics] = ...
    transfer_scatter_add_anisotropic(kVals, edges, E, ET, EF, HPOL, HDIR, HT, ...
                                     weight_k, centroidX_k, centroidY_k, ...
                                     mu1, mu3, nu, D, t, t0)
% transfer_scatter_add_anisotropic: Full anisotropic EDQNM transfer computation
%
% Computes 6 coupled transfer terms for anisotropic turbulence:
%   1. S_NL_ISO  - Isotropic energy transfer
%   2. S_NL_DIR  - Directional anisotropy transfer
%   3. S_NL_POL  - Poloidal anisotropy transfer
%   4. ST_NL_ISO - Scalar isotropic transfer
%   5. ST_NL_DIR - Scalar directional transfer
%   6. SF_NL     - Flux transfer
%
% All 6 computed in single k-slice loop for efficiency and parallelizability.
%
% K-SLICE INTEGRATION: Integrates each triad once via k-slice.
% Cyclic symmetry of the integration domain ensures correctness.
% Delta correction at each triad enforces local energy conservation.
%
% PARALLELIZATION STRATEGY (OpenMP/MPI):
%   - Parallelize outer loop over kj
%   - Each thread maintains 6 separate dE_local arrays with Kahan compensation
%   - Combine thread-local arrays at end using Kahan summation
%
% Inputs:
%   kVals, edges          - Wavenumber grid (kLength x 1)
%   E, ET, EF             - Energy spectra (kLength x 1)
%   HPOL, HDIR, HT        - Anisotropy coefficients (kLength x 1)
%   weight_k, centroids   - FV integration weights and centroids (kLength x kLength x kLength)
%   mu1, mu3              - Eddy damping coefficients (kLength x 1)
%   nu, D                 - Viscosity and diffusivity (scalars)
%   t, t0                 - Integration time and reference time (scalars)
%
% Outputs:
%   S_NL_ISO, S_NL_DIR, S_NL_POL - Energy transfer rates (kLength x 1)
%   ST_NL_ISO, ST_NL_DIR         - Scalar transfer rates (kLength x 1)
%   SF_NL                        - Flux transfer rates (kLength x 1)
%   diagnostics                  - Struct with diagnostics

kVals = double(kVals(:));
edges = double(edges(:));
E     = double(E(:));
ET    = double(ET(:));
EF    = double(EF(:));
HPOL  = double(HPOL(:));
HDIR  = double(HDIR(:));
HT    = double(HT(:));
mu1   = double(mu1(:));
mu3   = double(mu3(:));
nu    = double(nu);
D     = double(D);
t     = double(t);
t0    = double(t0);

weight_k    = double(weight_k);
centroidX_k = double(centroidX_k);
centroidY_k = double(centroidY_k);

kLength = length(kVals);
dk = diff(edges(:));
krep = sqrt(edges(1:end-1).*edges(2:end));  % FV reps

% Precompute spectral densities and interpolation data
E0  = E ./ (4*pi*kVals.^2);
E0s = max(E0, realmin);
E0T = ET ./ (4*pi*kVals.^2);
E0Ts = max(E0T, realmin);

logk  = log(kVals);
logE0 = log(E0s);
logE0T = log(E0Ts);
logEF = log(max(EF, realmin));
logHPOL = log(max(abs(HPOL), realmin));
logHDIR = log(max(abs(HDIR), realmin));
logHT = log(max(abs(HT), realmin));

% Sign arrays for anisotropy coefficients (preserve sign after log)
signHPOL = sign(HPOL);
signHDIR = sign(HDIR);
signHT = sign(HT);

% mu interpolation data
mu1s = max(mu1, realmin);
mu3s = max(mu3, realmin);
logmu1 = log(mu1s);
logmu3 = log(mu3s);

% Energy increment accumulators (6 separate arrays with Kahan compensation)
dE_iso = zeros(kLength,1);  cE_iso = zeros(kLength,1);
dE_dir = zeros(kLength,1);  cE_dir = zeros(kLength,1);
dE_pol = zeros(kLength,1);  cE_pol = zeros(kLength,1);
dET_iso = zeros(kLength,1); cET_iso = zeros(kLength,1);
dET_dir = zeros(kLength,1); cET_dir = zeros(kLength,1);
dEF = zeros(kLength,1);     cEF = zeros(kLength,1);

% Diagnostics
numTriadsUsed = 0;
maxTriadResidual_iso = 0.0;
maxTriadResidual_dir = 0.0;
maxTriadResidual_pol = 0.0;
maxTriadResidual_T_iso = 0.0;
maxTriadResidual_T_dir = 0.0;
maxTriadResidual_F = 0.0;

% ===========================================================================
% K-SLICE LOOP: Parallelizable over kj
% ===========================================================================
for kj = 1:kLength
    kstar = krep(kj);
    k2 = kstar^2;

    for pj = 1:kLength
        for qj = pj:kLength
            wk = weight_k(pj,qj,kj);
            if wk <= 0, continue; end

            dv_local = wk * dk(kj);

            % Evaluate kernel at (kstar, pstar, qstar) from k-slice centroids
            pstar = centroidX_k(pj,qj,kj);
            qstar = centroidY_k(pj,qj,kj);

            % Compute all 6 transfer increments for this triad
            [dEk_iso, dEp_iso, dEq_iso, ...
             dEk_dir, dEp_dir, dEq_dir, ...
             dEk_pol, dEp_pol, dEq_pol, ...
             dETk_iso, dETp_iso, dETq_iso, ...
             dETk_dir, dETp_dir, dETq_dir, ...
             dEFk, dEFp, dEFq] = ...
                compute_all_triad_increments(kstar, pstar, qstar, dv_local, kj, pj, qj, ...
                    logk, logE0, logE0T, logEF, ...
                    logHPOL, signHPOL, logHDIR, signHDIR, logHT, signHT, ...
                    logmu1, logmu3, nu, D, t, t0);

            % Scatter-add to bins using Kahan summation
            % S_NL_ISO
            [dE_iso(kj), cE_iso(kj)] = kahan_add(dE_iso(kj), cE_iso(kj), dEk_iso);
            [dE_iso(pj), cE_iso(pj)] = kahan_add(dE_iso(pj), cE_iso(pj), dEp_iso);
            [dE_iso(qj), cE_iso(qj)] = kahan_add(dE_iso(qj), cE_iso(qj), dEq_iso);

            % S_NL_DIR
            [dE_dir(kj), cE_dir(kj)] = kahan_add(dE_dir(kj), cE_dir(kj), dEk_dir);
            [dE_dir(pj), cE_dir(pj)] = kahan_add(dE_dir(pj), cE_dir(pj), dEp_dir);
            [dE_dir(qj), cE_dir(qj)] = kahan_add(dE_dir(qj), cE_dir(qj), dEq_dir);

            % S_NL_POL
            [dE_pol(kj), cE_pol(kj)] = kahan_add(dE_pol(kj), cE_pol(kj), dEk_pol);
            [dE_pol(pj), cE_pol(pj)] = kahan_add(dE_pol(pj), cE_pol(pj), dEp_pol);
            [dE_pol(qj), cE_pol(qj)] = kahan_add(dE_pol(qj), cE_pol(qj), dEq_pol);

            % ST_NL_ISO
            [dET_iso(kj), cET_iso(kj)] = kahan_add(dET_iso(kj), cET_iso(kj), dETk_iso);
            [dET_iso(pj), cET_iso(pj)] = kahan_add(dET_iso(pj), cET_iso(pj), dETp_iso);
            [dET_iso(qj), cET_iso(qj)] = kahan_add(dET_iso(qj), cET_iso(qj), dETq_iso);

            % ST_NL_DIR
            [dET_dir(kj), cET_dir(kj)] = kahan_add(dET_dir(kj), cET_dir(kj), dETk_dir);
            [dET_dir(pj), cET_dir(pj)] = kahan_add(dET_dir(pj), cET_dir(pj), dETp_dir);
            [dET_dir(qj), cET_dir(qj)] = kahan_add(dET_dir(qj), cET_dir(qj), dETq_dir);

            % SF_NL
            [dEF(kj), cEF(kj)] = kahan_add(dEF(kj), cEF(kj), dEFk);
            [dEF(pj), cEF(pj)] = kahan_add(dEF(pj), cEF(pj), dEFp);
            [dEF(qj), cEF(qj)] = kahan_add(dEF(qj), cEF(qj), dEFq);

            % Track max residuals
            maxTriadResidual_iso = max(maxTriadResidual_iso, abs(dEk_iso + dEp_iso + dEq_iso));
            maxTriadResidual_dir = max(maxTriadResidual_dir, abs(dEk_dir + dEp_dir + dEq_dir));
            maxTriadResidual_pol = max(maxTriadResidual_pol, abs(dEk_pol + dEp_pol + dEq_pol));
            maxTriadResidual_T_iso = max(maxTriadResidual_T_iso, abs(dETk_iso + dETp_iso + dETq_iso));
            maxTriadResidual_T_dir = max(maxTriadResidual_T_dir, abs(dETk_dir + dETp_dir + dETq_dir));
            maxTriadResidual_F = max(maxTriadResidual_F, abs(dEFk + dEFp + dEFq));

            numTriadsUsed = numTriadsUsed + 1;
        end
    end
end

% Convert energy increments to transfer rates
S_NL_ISO = dE_iso ./ dk;
S_NL_DIR = dE_dir ./ dk;
S_NL_POL = dE_pol ./ dk;
ST_NL_ISO = dET_iso ./ dk;
ST_NL_DIR = dET_dir ./ dk;
SF_NL = dEF ./ dk;

% Package diagnostics
diagnostics.numTriadsUsed = numTriadsUsed;
diagnostics.FV_total_transfer_iso = sum(dE_iso);
diagnostics.FV_total_transfer_dir = sum(dE_dir);
diagnostics.FV_total_transfer_pol = sum(dE_pol);
diagnostics.FV_total_transfer_T_iso = sum(dET_iso);
diagnostics.FV_total_transfer_T_dir = sum(dET_dir);
diagnostics.FV_total_transfer_F = sum(dEF);
diagnostics.maxTriadResidual_iso = maxTriadResidual_iso;
diagnostics.maxTriadResidual_dir = maxTriadResidual_dir;
diagnostics.maxTriadResidual_pol = maxTriadResidual_pol;
diagnostics.maxTriadResidual_T_iso = maxTriadResidual_T_iso;
diagnostics.maxTriadResidual_T_dir = maxTriadResidual_T_dir;
diagnostics.maxTriadResidual_F = maxTriadResidual_F;

end

% ============================================================================
% MAIN COMPUTATION: All 6 triad increments
% ============================================================================
function [dEk_iso, dEp_iso, dEq_iso, ...
          dEk_dir, dEp_dir, dEq_dir, ...
          dEk_pol, dEp_pol, dEq_pol, ...
          dETk_iso, dETp_iso, dETq_iso, ...
          dETk_dir, dETp_dir, dETq_dir, ...
          dEFk, dEFp, dEFq] = ...
    compute_all_triad_increments(kstar, pstar, qstar, dv, kj, pj, qj, ...
        logk, logE0, logE0T, logEF, ...
        logHPOL, signHPOL, logHDIR, signHDIR, logHT, signHT, ...
        logmu1, logmu3, nu, D, t, t0)

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

% Triad geometry: SAME for all legs (fixed k, p, q configuration)
k = kstar; p = pstar; q = qstar;
k2 = k^2; p2 = p^2; q2 = q^2;
x = (k2 + p2 - q2) / (2*k*p);  % cos(angle between k and p)
y = (k2 + q2 - p2) / (2*k*q);  % cos(angle between k and q)
z = (p2 + q2 - k2) / (2*p*q);  % cos(angle between p and q)

% Interpolate all spectral quantities at triad centroids
E0k = interp_log(k, logk, logE0);
E0p = interp_log(p, logk, logE0);
E0q = interp_log(q, logk, logE0);

E0Tk = interp_log(k, logk, logE0T);
E0Tp = interp_log(p, logk, logE0T);
E0Tq = interp_log(q, logk, logE0T);

% Reconstruct E, ET from E0, E0T
Ek = 4*pi*k2*E0k;
Ep = 4*pi*p2*E0p;
Eq = 4*pi*q2*E0q;
ETk = 4*pi*k2*E0Tk;
ETp = 4*pi*p2*E0Tp;
ETq = 4*pi*q2*E0Tq;

EFk = interp_log(k, logk, logEF);
EFp = interp_log(p, logk, logEF);
EFq = interp_log(q, logk, logEF);

HPOLk = interp_log_signed(k, logk, logHPOL, signHPOL);
HPOLp = interp_log_signed(p, logk, logHPOL, signHPOL);
HPOLq = interp_log_signed(q, logk, logHPOL, signHPOL);

HDIRk = interp_log_signed(k, logk, logHDIR, signHDIR);
HDIRp = interp_log_signed(p, logk, logHDIR, signHDIR);
HDIRq = interp_log_signed(q, logk, logHDIR, signHDIR);

HTk = interp_log_signed(k, logk, logHT, signHT);
HTp = interp_log_signed(p, logk, logHT, signHT);
HTq = interp_log_signed(q, logk, logHT, signHT);

mu1_k = interp_log(k, logk, logmu1);
mu1_p = interp_log(p, logk, logmu1);
mu1_q = interp_log(q, logk, logmu1);

mu3_k = interp_log(k, logk, logmu3);
mu3_p = interp_log(p, logk, logmu3);
mu3_q = interp_log(q, logk, logmu3);

% Compute theta functions
theta_val = theta(nu, k, p, q, kj, pj, qj, t, mu1_k, mu1_p, mu1_q);
thetaT_val = thetaT(nu, D, k, p, q, t-t0, mu3_q);
thetaF_kpq = thetaF(nu, D, k, p, q, t-t0, mu3_p, mu3_q);
thetaF_pkq = thetaF(nu, D, p, k, q, t-t0, mu3_k, mu3_q);

% ============================================================================
% 1. S_NL_ISO (Isotropic energy transfer)
% ============================================================================
[dEk_iso, dEp_iso, dEq_iso] = compute_S_NL_ISO(...
    E0k, E0p, E0q, k, p, q, x, y, z, theta_val, dv);

% ============================================================================
% 2. S_NL_DIR (Directional anisotropy transfer)
% ============================================================================
[dEk_dir, dEp_dir, dEq_dir] = compute_S_NL_DIR(...
    E0k, E0p, E0q, HPOLk, HPOLp, HPOLq, HDIRk, HDIRp, HDIRq, ...
    k, p, q, x, y, z, theta_val, dv);

% ============================================================================
% 3. S_NL_POL (Poloidal anisotropy transfer)
% ============================================================================
[dEk_pol, dEp_pol, dEq_pol] = compute_S_NL_POL(...
    E0k, E0p, E0q, HPOLk, HPOLp, HPOLq, HDIRk, HDIRp, HDIRq, ...
    k, p, q, x, y, z, theta_val, dv);

% ============================================================================
% 4. ST_NL_ISO (Scalar isotropic transfer)
% ============================================================================
[dETk_iso, dETp_iso, dETq_iso] = compute_ST_NL_ISO(...
    Ek, Ep, Eq, k, p, q, E0Tk, E0Tp, E0Tq, k2, p2, y, thetaT_val, dv);

% ============================================================================
% 5. ST_NL_DIR (Scalar directional transfer)
% ============================================================================
[dETk_dir, dETp_dir, dETq_dir] = compute_ST_NL_DIR(...
    E0k, E0p, E0q, E0Tk, E0Tp, E0Tq, ...
    HPOLk, HPOLp, HPOLq, HDIRk, HDIRp, HDIRq, HTk, HTp, HTq, ...
    k, p, q, k2, p2, q2, x, y, z, thetaT_val, dv);

% ============================================================================
% 6. SF_NL (Flux transfer)
% ============================================================================
[dEFk, dEFp, dEFq] = compute_SF_NL(...
    EFk, EFp, EFq, E0k, E0p, E0q, ...
    k, p, q, k2, x, y, z, thetaF_kpq, thetaF_pkq, dv);

end

% ============================================================================
% 1. S_NL_ISO
% ============================================================================
function [dEk, dEp, dEq] = compute_S_NL_ISO(E0k, E0p, E0q, k, p, q, x, y, z, theta_val, dv)
% S_NL_ISO: theta * 16*pi^2 * p^2*k^2*q * (xy+z^3) * kernel1b
% kernel1b = E0(pj)*(E0(qj)-E0(kj)) with cyclic permutations

geom = theta_val * 16.0 * pi^2 * p^2 * k^2 * q * (x*y + z^3);

% k-leg: 0.5*(E0q*(E0p-E0k) + E0p*(E0q-E0k))
Sk_raw = geom * 0.5 * (E0q*(E0p-E0k) + E0p*(E0q-E0k));

% p-leg: 0.5*(E0k*(E0q-E0p) + E0q*(E0k-E0p))
Sp_raw = geom * 0.5 * (E0k*(E0q-E0p) + E0q*(E0k-E0p));

% q-leg: 0.5*(E0k*(E0p-E0q) + E0p*(E0k-E0q))
Sq_raw = geom * 0.5 * (E0k*(E0p-E0q) + E0p*(E0k-E0q));

% Delta correction and exact zero projection
delta = (Sk_raw + Sp_raw + Sq_raw) / 3.0;
dEk = (Sk_raw - delta) * dv;
dEp = (Sp_raw - delta) * dv;
dEq = (Sq_raw - delta) * dv;
s = dEk + dEp + dEq;
dEk = dEk - s/3.0;
dEp = dEp - s/3.0;
dEq = dEq - s/3.0;
end

% ============================================================================
% 2. S_NL_DIR
% ============================================================================
function [dEk, dEp, dEq] = compute_S_NL_DIR(...
    E0k, E0p, E0q, HPOLk, HPOLp, HPOLq, HDIRk, HDIRp, HDIRq, ...
    k, p, q, x, y, z, theta_val, dv)
% S_NL_DIR with kernels 21-25

k2 = k^2; p2 = p^2; q2 = q^2;
y2 = y^2; z2 = z^2;
xy_z3 = x*y + z^3;

% k-leg kernels (original spectral ordering)
kernel21_k = 0.5 * (E0q*(E0p-E0k)*HPOLq + E0p*(E0q-E0k)*HPOLp);
kernel22_k = 0.5 * (E0q*E0p*HPOLp + E0p*E0q*HPOLq);
kernel23_k = 0.5 * (E0q*(E0p-E0k)*HDIRq + E0p*(E0q-E0k)*HDIRp);
kernel24_k = 0.5 * (E0q*E0p*HDIRp + E0p*E0q*HDIRq);
kernel25_k = 0.5 * (E0q*E0k*HDIRk + E0p*E0k*HDIRk);

Sk_raw = theta_val * 4.0 * pi^2 * p2 * k2 * q * (...
    (y2 - 1.0) * xy_z3 * kernel21_k + ...
    z * (1.0 - z2)^2 * kernel22_k) + ...
    theta_val * 8.0 * pi^2 * p2 * k2 * q * xy_z3 * (...
    (3.0*y2 - 1.0) * kernel23_k + ...
    (3.0*z2 - 1.0) * kernel24_k - ...
    2.0 * kernel25_k);

% p-leg kernels (cyclic permutation: k→p, p→q, q→k)
kernel21_p = 0.5 * (E0k*(E0q-E0p)*HPOLk + E0q*(E0k-E0p)*HPOLq);
kernel22_p = 0.5 * (E0k*E0q*HPOLq + E0q*E0k*HPOLk);
kernel23_p = 0.5 * (E0k*(E0q-E0p)*HDIRk + E0q*(E0k-E0p)*HDIRq);
kernel24_p = 0.5 * (E0k*E0q*HDIRq + E0q*E0k*HDIRk);
kernel25_p = 0.5 * (E0k*E0p*HDIRp + E0q*E0p*HDIRp);

Sp_raw = theta_val * 4.0 * pi^2 * p2 * k2 * q * (...
    (y2 - 1.0) * xy_z3 * kernel21_p + ...
    z * (1.0 - z2)^2 * kernel22_p) + ...
    theta_val * 8.0 * pi^2 * p2 * k2 * q * xy_z3 * (...
    (3.0*y2 - 1.0) * kernel23_p + ...
    (3.0*z2 - 1.0) * kernel24_p - ...
    2.0 * kernel25_p);

% q-leg kernels (cyclic permutation: k→q, p→k, q→p)
kernel21_q = 0.5 * (E0p*(E0k-E0q)*HPOLp + E0k*(E0p-E0q)*HPOLk);
kernel22_q = 0.5 * (E0p*E0k*HPOLk + E0k*E0p*HPOLp);
kernel23_q = 0.5 * (E0p*(E0k-E0q)*HDIRp + E0k*(E0p-E0q)*HDIRk);
kernel24_q = 0.5 * (E0p*E0k*HDIRk + E0k*E0p*HDIRp);
kernel25_q = 0.5 * (E0p*E0q*HDIRq + E0k*E0q*HDIRq);

Sq_raw = theta_val * 4.0 * pi^2 * p2 * k2 * q * (...
    (y2 - 1.0) * xy_z3 * kernel21_q + ...
    z * (1.0 - z2)^2 * kernel22_q) + ...
    theta_val * 8.0 * pi^2 * p2 * k2 * q * xy_z3 * (...
    (3.0*y2 - 1.0) * kernel23_q + ...
    (3.0*z2 - 1.0) * kernel24_q - ...
    2.0 * kernel25_q);

% Delta correction and exact zero projection
delta = (Sk_raw + Sp_raw + Sq_raw) / 3.0;
dEk = (Sk_raw - delta) * dv;
dEp = (Sp_raw - delta) * dv;
dEq = (Sq_raw - delta) * dv;
s = dEk + dEp + dEq;
dEk = dEk - s/3.0;
dEp = dEp - s/3.0;
dEq = dEq - s/3.0;
end

% ============================================================================
% 3. S_NL_POL
% ============================================================================
function [dEk, dEp, dEq] = compute_S_NL_POL(...
    E0k, E0p, E0q, HPOLk, HPOLp, HPOLq, HDIRk, HDIRp, HDIRq, ...
    k, p, q, x, y, z, theta_val, dv)
% S_NL_POL with kernels 31-37

k2 = k^2; p2 = p^2;
y2 = y^2; z2 = z^2;
xy_z3 = x*y + z^3;

% k-leg kernels
kernel31_k = 0.5 * (E0q*E0p*HPOLp + E0p*E0q*HPOLq);
kernel32_k = 0.5 * (E0q*E0k*HPOLk + E0p*E0k*HPOLk);
kernel33_k = 0.5 * (E0q*(E0p-E0k)*HPOLq + E0p*(E0q-E0k)*HPOLp);
kernel34_k = kernel31_k;  % Same as kernel31
kernel35_k = 0.5 * (E0q*E0k*HPOLq + E0p*E0k*HPOLp);
kernel36_k = 0.5 * (E0q*(E0p-E0k)*HDIRq + E0p*(E0q-E0k)*HDIRp);
kernel37_k = 0.5 * (E0q*E0p*HDIRp + E0p*E0q*HDIRq);

Sk_raw = theta_val * 4.0 * pi^2 * p2 * k2 * q * (...
    xy_z3 * ((1.0 + z2) * kernel31_k - 4.0 * kernel32_k) + ...
    z * (z2 - 1.0) * (1.0 + y2) * kernel33_k + ...
    2.0 * z * (z2 - y2) * kernel34_k + ...
    2.0 * y * x * (z2 - 1.0) * kernel35_k) + ...
    theta_val * 24.0 * pi^2 * p2 * k2 * q * z * (z2 - 1.0) * (...
    (y2 - 1.0) * kernel36_k + (z2 - 1.0) * kernel37_k);

% p-leg kernels (cyclic permutation)
kernel31_p = 0.5 * (E0k*E0q*HPOLq + E0q*E0k*HPOLk);
kernel32_p = 0.5 * (E0k*E0p*HPOLp + E0q*E0p*HPOLp);
kernel33_p = 0.5 * (E0k*(E0q-E0p)*HPOLk + E0q*(E0k-E0p)*HPOLq);
kernel34_p = kernel31_p;
kernel35_p = 0.5 * (E0k*E0p*HPOLk + E0q*E0p*HPOLq);
kernel36_p = 0.5 * (E0k*(E0q-E0p)*HDIRk + E0q*(E0k-E0p)*HDIRq);
kernel37_p = 0.5 * (E0k*E0q*HDIRq + E0q*E0k*HDIRk);

Sp_raw = theta_val * 4.0 * pi^2 * p2 * k2 * q * (...
    xy_z3 * ((1.0 + z2) * kernel31_p - 4.0 * kernel32_p) + ...
    z * (z2 - 1.0) * (1.0 + y2) * kernel33_p + ...
    2.0 * z * (z2 - y2) * kernel34_p + ...
    2.0 * y * x * (z2 - 1.0) * kernel35_p) + ...
    theta_val * 24.0 * pi^2 * p2 * k2 * q * z * (z2 - 1.0) * (...
    (y2 - 1.0) * kernel36_p + (z2 - 1.0) * kernel37_p);

% q-leg kernels (cyclic permutation)
kernel31_q = 0.5 * (E0p*E0k*HPOLk + E0k*E0p*HPOLp);
kernel32_q = 0.5 * (E0p*E0q*HPOLq + E0k*E0q*HPOLq);
kernel33_q = 0.5 * (E0p*(E0k-E0q)*HPOLp + E0k*(E0p-E0q)*HPOLk);
kernel34_q = kernel31_q;
kernel35_q = 0.5 * (E0p*E0q*HPOLp + E0k*E0q*HPOLk);
kernel36_q = 0.5 * (E0p*(E0k-E0q)*HDIRp + E0k*(E0p-E0q)*HDIRk);
kernel37_q = 0.5 * (E0p*E0k*HDIRk + E0k*E0p*HDIRp);

Sq_raw = theta_val * 4.0 * pi^2 * p2 * k2 * q * (...
    xy_z3 * ((1.0 + z2) * kernel31_q - 4.0 * kernel32_q) + ...
    z * (z2 - 1.0) * (1.0 + y2) * kernel33_q + ...
    2.0 * z * (z2 - y2) * kernel34_q + ...
    2.0 * y * x * (z2 - 1.0) * kernel35_q) + ...
    theta_val * 24.0 * pi^2 * p2 * k2 * q * z * (z2 - 1.0) * (...
    (y2 - 1.0) * kernel36_q + (z2 - 1.0) * kernel37_q);

% Delta correction and exact zero projection
delta = (Sk_raw + Sp_raw + Sq_raw) / 3.0;
dEk = (Sk_raw - delta) * dv;
dEp = (Sp_raw - delta) * dv;
dEq = (Sq_raw - delta) * dv;
s = dEk + dEp + dEq;
dEk = dEk - s/3.0;
dEp = dEp - s/3.0;
dEq = dEq - s/3.0;
end

% ============================================================================
% 4. ST_NL_ISO
% ============================================================================
function [dETk, dETp, dETq] = compute_ST_NL_ISO(...
    Ek, Ep, Eq, k, p, q, E0Tk, E0Tp, E0Tq, k2, p2, y, thetaT_val, dv)
% ST_NL_ISO: thetaT * k/p/q * (1-y^2) * kernel4
% kernel4 = Eq * (k^2*E0T(p) - p^2*E0T(k))

geom = thetaT_val * k / (p*q) * (1.0 - y^2);
q2 = q^2;

% k-leg: Eq * (k^2*E0Tp - p^2*E0Tk)
STk_raw = geom * Eq * (k2*E0Tp - p2*E0Tk);

% p-leg (cyclic): Ek * (p^2*E0Tq - q^2*E0Tp)
STp_raw = geom * Ek * (p2*E0Tq - q2*E0Tp);

% q-leg (cyclic): Ep * (q^2*E0Tk - k^2*E0Tq)
STq_raw = geom * Ep * (q2*E0Tk - k2*E0Tq);

% Delta correction and exact zero projection
delta = (STk_raw + STp_raw + STq_raw) / 3.0;
dETk = (STk_raw - delta) * dv;
dETp = (STp_raw - delta) * dv;
dETq = (STq_raw - delta) * dv;
s = dETk + dETp + dETq;
dETk = dETk - s/3.0;
dETp = dETp - s/3.0;
dETq = dETq - s/3.0;
end

% ============================================================================
% 5. ST_NL_DIR
% ============================================================================
function [dETk, dETp, dETq] = compute_ST_NL_DIR(...
    E0k, E0p, E0q, E0Tk, E0Tp, E0Tq, ...
    HPOLk, HPOLp, HPOLq, HDIRk, HDIRp, HDIRq, HTk, HTp, HTq, ...
    k, p, q, k2, p2, q2, x, y, z, thetaT_val, dv)
% ST_NL_DIR with kernels 51-54

y2 = y^2; z2 = z^2;
xy_z = x*y + z;

% k-leg kernels
kernel51_k = 0.5 * (E0q*(E0Tp-E0Tk)*HPOLq + E0p*(E0Tq-E0Tk)*HPOLp);
kernel52_k = 0.5 * (E0q*(E0Tp-E0Tk)*HDIRq + E0p*(E0Tq-E0Tk)*HDIRp);
kernel53_k = 0.5 * (E0q*E0Tp*HTp + E0p*E0Tq*HTq);
kernel54_k = 0.5 * (E0q*2.0*E0Tk*HTk + E0p*2.0*E0Tk*HTk);

STk_raw = 4.0 * thetaT_val * pi^2 * k2 * p2 * q * xy_z * (y2 - 1.0) * kernel51_k + ...
    8.0 * thetaT_val * pi^2 * k2 * p2 * q * xy_z * (3.0*y2 - 1.0) * kernel52_k + ...
    8.0 * thetaT_val * pi^2 * k2 * p2 * q * xy_z * ((3.0*z2 - 1.0) * kernel53_k - kernel54_k);

% p-leg kernels (cyclic permutation: k→p, p→q, q→k)
kernel51_p = 0.5 * (E0k*(E0Tq-E0Tp)*HPOLk + E0q*(E0Tk-E0Tp)*HPOLq);
kernel52_p = 0.5 * (E0k*(E0Tq-E0Tp)*HDIRk + E0q*(E0Tk-E0Tp)*HDIRq);
kernel53_p = 0.5 * (E0k*E0Tq*HTq + E0q*E0Tk*HTk);
kernel54_p = 0.5 * (E0k*2.0*E0Tp*HTp + E0q*2.0*E0Tp*HTp);

STp_raw = 4.0 * thetaT_val * pi^2 * k2 * p2 * q * xy_z * (y2 - 1.0) * kernel51_p + ...
    8.0 * thetaT_val * pi^2 * k2 * p2 * q * xy_z * (3.0*y2 - 1.0) * kernel52_p + ...
    8.0 * thetaT_val * pi^2 * k2 * p2 * q * xy_z * ((3.0*z2 - 1.0) * kernel53_p - kernel54_p);

% q-leg kernels (cyclic permutation: k→q, p→k, q→p)
kernel51_q = 0.5 * (E0p*(E0Tk-E0Tq)*HPOLp + E0k*(E0Tp-E0Tq)*HPOLk);
kernel52_q = 0.5 * (E0p*(E0Tk-E0Tq)*HDIRp + E0k*(E0Tp-E0Tq)*HDIRk);
kernel53_q = 0.5 * (E0p*E0Tk*HTk + E0k*E0Tp*HTp);
kernel54_q = 0.5 * (E0p*2.0*E0Tq*HTq + E0k*2.0*E0Tq*HTq);

STq_raw = 4.0 * thetaT_val * pi^2 * k2 * p2 * q * xy_z * (y2 - 1.0) * kernel51_q + ...
    8.0 * thetaT_val * pi^2 * k2 * p2 * q * xy_z * (3.0*y2 - 1.0) * kernel52_q + ...
    8.0 * thetaT_val * pi^2 * k2 * p2 * q * xy_z * ((3.0*z2 - 1.0) * kernel53_q - kernel54_q);

% Delta correction and exact zero projection
delta = (STk_raw + STp_raw + STq_raw) / 3.0;
dETk = (STk_raw - delta) * dv;
dETp = (STp_raw - delta) * dv;
dETq = (STq_raw - delta) * dv;
s = dETk + dETp + dETq;
dETk = dETk - s/3.0;
dETp = dETp - s/3.0;
dETq = dETq - s/3.0;
end

% ============================================================================
% 6. SF_NL
% ============================================================================
function [dEFk, dEFp, dEFq] = compute_SF_NL(...
    EFk, EFp, EFq, E0k, E0p, E0q, ...
    k, p, q, k2, x, y, z, thetaF_kpq, thetaF_pkq, dv)
% SF_NL with kernels 61-66

y2 = y^2; z2 = z^2; y3 = y2*y;

% k-leg kernels
kernel61_k = E0p*EFq;
kernel62_k = E0p*EFk;
kernel63_k = E0k*EFp;
kernel64_k = E0k*EFq;
kernel65_k = E0q*EFp;
kernel66_k = E0q*EFk;

SFk_raw = 4.0 * pi^2 * thetaF_kpq * k2 * p * q * (...
    k * kernel61_k * (1.0 + y2 - z2 - x*y*z - 2.0*y2*z2) - ...
    2.0 * q * (y3 + x*z) * kernel62_k) + ...
    4.0 * pi^2 * thetaF_pkq * k2 * p * q * (...
    q * z * (2.0*x*y2 + y*z - x) * kernel63_k - ...
    p * y * (x + y*z) * kernel64_k + ...
    k * ((1.0 - y2 + z2 - x*y*z - 2.0*y2*z2) * kernel65_k - 2.0*(1.0 - y2) * kernel66_k));

% p-leg kernels (cyclic permutation: k→p, p→q, q→k)
kernel61_p = E0q*EFk;
kernel62_p = E0q*EFp;
kernel63_p = E0p*EFq;
kernel64_p = E0p*EFk;
kernel65_p = E0k*EFq;
kernel66_p = E0k*EFp;

SFp_raw = 4.0 * pi^2 * thetaF_kpq * k2 * p * q * (...
    p * kernel61_p * (1.0 + y2 - z2 - x*y*z - 2.0*y2*z2) - ...
    2.0 * k * (y3 + x*z) * kernel62_p) + ...
    4.0 * pi^2 * thetaF_pkq * k2 * p * q * (...
    k * z * (2.0*x*y2 + y*z - x) * kernel63_p - ...
    q * y * (x + y*z) * kernel64_p + ...
    p * ((1.0 - y2 + z2 - x*y*z - 2.0*y2*z2) * kernel65_p - 2.0*(1.0 - y2) * kernel66_p));

% q-leg kernels (cyclic permutation: k→q, p→k, q→p)
kernel61_q = E0k*EFp;
kernel62_q = E0k*EFq;
kernel63_q = E0q*EFk;
kernel64_q = E0q*EFp;
kernel65_q = E0p*EFk;
kernel66_q = E0p*EFq;

SFq_raw = 4.0 * pi^2 * thetaF_kpq * k2 * p * q * (...
    q * kernel61_q * (1.0 + y2 - z2 - x*y*z - 2.0*y2*z2) - ...
    2.0 * p * (y3 + x*z) * kernel62_q) + ...
    4.0 * pi^2 * thetaF_pkq * k2 * p * q * (...
    p * z * (2.0*x*y2 + y*z - x) * kernel63_q - ...
    k * y * (x + y*z) * kernel64_q + ...
    q * ((1.0 - y2 + z2 - x*y*z - 2.0*y2*z2) * kernel65_q - 2.0*(1.0 - y2) * kernel66_q));

% Delta correction and exact zero projection
delta = (SFk_raw + SFp_raw + SFq_raw) / 3.0;
dEFk = (SFk_raw - delta) * dv;
dEFp = (SFp_raw - delta) * dv;
dEFq = (SFq_raw - delta) * dv;
s = dEFk + dEFp + dEFq;
dEFk = dEFk - s/3.0;
dEFp = dEFp - s/3.0;
dEFq = dEFq - s/3.0;
end

% ============================================================================
% INTERPOLATION HELPER FUNCTIONS
% ============================================================================

function val = interp_log(k_eval, logk, logval)
% Log-log linear interpolation
logk_eval = log(k_eval);
logval_interp = interp1(logk, logval, logk_eval, 'linear', 'extrap');
val = exp(logval_interp);
end

function val = interp_log_signed(k_eval, logk, logval, signval)
% Log-log interpolation preserving sign
logk_eval = log(k_eval);
sign_interp = interp1(logk, signval, logk_eval, 'linear', 'extrap');
logval_interp = interp1(logk, logval, logk_eval, 'linear', 'extrap');
val = sign_interp * exp(logval_interp);
end

% ============================================================================
% KAHAN SUMMATION
% ============================================================================

function [sum_new, c_new] = kahan_add(sum_old, c_old, x)
y = x - c_old;
t = sum_old + y;
c_new = (t - sum_old) - y;
sum_new = t;
end
