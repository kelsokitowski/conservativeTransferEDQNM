function [S_NL_ISO, S_NL_DIR, S_NL_POL, ST_NL_ISO, ST_NL_DIR, SF_NL, diagnostics] = ...
    transfer_scatter_add_anisotropic(kVals, edges, E, ET, EF, HPOL, HDIR, HT, ...
                                     weight_k, centroidX_k, centroidY_k, ...
                                     mu1, mu3, nu, D, t, t0)
% transfer_scatter_add_anisotropic: Full anisotropic EDQNM transfer computation
%
% Computes 6 coupled transfer terms for anisotropic turbulence:
%   1. S_NL_ISO  - Isotropic energy transfer (kernel1)
%   2. S_NL_DIR  - Directional anisotropy transfer (kernel2)
%   3. S_NL_POL  - Poloidal anisotropy transfer (kernel3)
%   4. ST_NL_ISO - Scalar isotropic transfer (kernel4)
%   5. ST_NL_DIR - Scalar directional transfer (kernel5)
%   6. SF_NL     - Flux transfer (kernel6)
%
% All 6 computed in single k-slice loop for efficiency and parallelizability.
% Each kernel function encapsulates full geometric expression.
% Symmetrization: Sk_raw = 0.5*(kernel(...term1...) + kernel(...term2...))

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
krep = sqrt(edges(1:end-1).*edges(2:end));

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

signHPOL = sign(HPOL);
signHDIR = sign(HDIR);
signHT = sign(HT);

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

    for pj = 1:kLength
        for qj = pj:kLength
            wk = weight_k(pj,qj,kj);
            if wk <= 0, continue; end

            dv_local = wk * dk(kj);

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

            % Scatter-add using Kahan summation
            [dE_iso(kj), cE_iso(kj)] = kahan_add(dE_iso(kj), cE_iso(kj), dEk_iso);
            [dE_iso(pj), cE_iso(pj)] = kahan_add(dE_iso(pj), cE_iso(pj), dEp_iso);
            [dE_iso(qj), cE_iso(qj)] = kahan_add(dE_iso(qj), cE_iso(qj), dEq_iso);

            [dE_dir(kj), cE_dir(kj)] = kahan_add(dE_dir(kj), cE_dir(kj), dEk_dir);
            [dE_dir(pj), cE_dir(pj)] = kahan_add(dE_dir(pj), cE_dir(pj), dEp_dir);
            [dE_dir(qj), cE_dir(qj)] = kahan_add(dE_dir(qj), cE_dir(qj), dEq_dir);

            [dE_pol(kj), cE_pol(kj)] = kahan_add(dE_pol(kj), cE_pol(kj), dEk_pol);
            [dE_pol(pj), cE_pol(pj)] = kahan_add(dE_pol(pj), cE_pol(pj), dEp_pol);
            [dE_pol(qj), cE_pol(qj)] = kahan_add(dE_pol(qj), cE_pol(qj), dEq_pol);

            [dET_iso(kj), cET_iso(kj)] = kahan_add(dET_iso(kj), cET_iso(kj), dETk_iso);
            [dET_iso(pj), cET_iso(pj)] = kahan_add(dET_iso(pj), cET_iso(pj), dETp_iso);
            [dET_iso(qj), cET_iso(qj)] = kahan_add(dET_iso(qj), cET_iso(qj), dETq_iso);

            [dET_dir(kj), cET_dir(kj)] = kahan_add(dET_dir(kj), cET_dir(kj), dETk_dir);
            [dET_dir(pj), cET_dir(pj)] = kahan_add(dET_dir(pj), cET_dir(pj), dETp_dir);
            [dET_dir(qj), cET_dir(qj)] = kahan_add(dET_dir(qj), cET_dir(qj), dETq_dir);

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

% Convert to transfer rates
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

k = kstar; p = pstar; q = qstar;
k2 = k^2; p2 = p^2; q2 = q^2;

% Interpolate all spectral quantities
E0k = interp_log(k, logk, logE0);
E0p = interp_log(p, logk, logE0);
E0q = interp_log(q, logk, logE0);

E0Tk = interp_log(k, logk, logE0T);
E0Tp = interp_log(p, logk, logE0T);
E0Tq = interp_log(q, logk, logE0T);

Ek = 4*pi*k2*E0k;
Ep = 4*pi*p2*E0p;
Eq = 4*pi*q2*E0q;

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

% ============================================================================
% 1. S_NL_ISO: kernel1
% ============================================================================
% ALL arguments fully permuted for each leg
% k-leg: Sk_raw = 0.5*(kernel1(k,p,q) + kernel1(k,q,p))
term1_k = kernel1(E0p, E0q, E0k, k, p, q, mu1_k, mu1_p, mu1_q, kj, pj, qj, nu, t);
term2_k = kernel1(E0q, E0p, E0k, k, q, p, mu1_k, mu1_q, mu1_p, kj, qj, pj, nu, t);
Sk_raw = 0.5 * (term1_k + term2_k);

% p-leg: Sp_raw = 0.5*(kernel1(p,q,k) + kernel1(p,k,q))
term1_p = kernel1(E0q, E0k, E0p, p, q, k, mu1_p, mu1_q, mu1_k, pj, qj, kj, nu, t);
term2_p = kernel1(E0k, E0q, E0p, p, k, q, mu1_p, mu1_k, mu1_q, pj, kj, qj, nu, t);
Sp_raw = 0.5 * (term1_p + term2_p);

% q-leg: Sq_raw = 0.5*(kernel1(q,k,p) + kernel1(q,p,k))
term1_q = kernel1(E0k, E0p, E0q, q, k, p, mu1_q, mu1_k, mu1_p, qj, kj, pj, nu, t);
term2_q = kernel1(E0p, E0k, E0q, q, p, k, mu1_q, mu1_p, mu1_k, qj, pj, kj, nu, t);
Sq_raw = 0.5 * (term1_q + term2_q);

delta = (Sk_raw + Sp_raw + Sq_raw) / 3.0;
dEk_iso = (Sk_raw - delta) * dv;
dEp_iso = (Sp_raw - delta) * dv;
dEq_iso = (Sq_raw - delta) * dv;
s = dEk_iso + dEp_iso + dEq_iso;
dEk_iso = dEk_iso - s/3.0;
dEp_iso = dEp_iso - s/3.0;
dEq_iso = dEq_iso - s/3.0;

% ============================================================================
% 2. S_NL_DIR: kernel2
% ============================================================================
% k-leg: Sk_raw = 0.5*(kernel2(k,p,q) + kernel2(k,q,p))
term1_k = kernel2(E0p, E0q, E0k, HPOLp, HPOLq, HDIRp, HDIRq, HDIRk, k, p, q, mu1_k, mu1_p, mu1_q, kj, pj, qj, nu, t);
term2_k = kernel2(E0q, E0p, E0k, HPOLq, HPOLp, HDIRq, HDIRp, HDIRk, k, q, p, mu1_k, mu1_q, mu1_p, kj, qj, pj, nu, t);
Sk_raw = 0.5 * (term1_k + term2_k);

% p-leg: Sp_raw = 0.5*(kernel2(p,q,k) + kernel2(p,k,q))
term1_p = kernel2(E0q, E0k, E0p, HPOLq, HPOLk, HDIRq, HDIRk, HDIRp, p, q, k, mu1_p, mu1_q, mu1_k, pj, qj, kj, nu, t);
term2_p = kernel2(E0k, E0q, E0p, HPOLk, HPOLq, HDIRk, HDIRq, HDIRp, p, k, q, mu1_p, mu1_k, mu1_q, pj, kj, qj, nu, t);
Sp_raw = 0.5 * (term1_p + term2_p);

% q-leg: Sq_raw = 0.5*(kernel2(q,k,p) + kernel2(q,p,k))
term1_q = kernel2(E0k, E0p, E0q, HPOLk, HPOLp, HDIRk, HDIRp, HDIRq, q, k, p, mu1_q, mu1_k, mu1_p, qj, kj, pj, nu, t);
term2_q = kernel2(E0p, E0k, E0q, HPOLp, HPOLk, HDIRp, HDIRk, HDIRq, q, p, k, mu1_q, mu1_p, mu1_k, qj, pj, kj, nu, t);
Sq_raw = 0.5 * (term1_q + term2_q);

delta = (Sk_raw + Sp_raw + Sq_raw) / 3.0;
dEk_dir = (Sk_raw - delta) * dv;
dEp_dir = (Sp_raw - delta) * dv;
dEq_dir = (Sq_raw - delta) * dv;
s = dEk_dir + dEp_dir + dEq_dir;
dEk_dir = dEk_dir - s/3.0;
dEp_dir = dEp_dir - s/3.0;
dEq_dir = dEq_dir - s/3.0;

% ============================================================================
% 3. S_NL_POL: kernel3
% ============================================================================
% k-leg: Sk_raw = 0.5*(kernel3(k,p,q) + kernel3(k,q,p))
term1_k = kernel3(E0p, E0q, E0k, HPOLp, HPOLq, HPOLk, HDIRp, HDIRq, k, p, q, mu1_k, mu1_p, mu1_q, kj, pj, qj, nu, t);
term2_k = kernel3(E0q, E0p, E0k, HPOLq, HPOLp, HPOLk, HDIRq, HDIRp, k, q, p, mu1_k, mu1_q, mu1_p, kj, qj, pj, nu, t);
Sk_raw = 0.5 * (term1_k + term2_k);

% p-leg: Sp_raw = 0.5*(kernel3(p,q,k) + kernel3(p,k,q))
term1_p = kernel3(E0q, E0k, E0p, HPOLq, HPOLk, HPOLp, HDIRq, HDIRk, p, q, k, mu1_p, mu1_q, mu1_k, pj, qj, kj, nu, t);
term2_p = kernel3(E0k, E0q, E0p, HPOLk, HPOLq, HPOLp, HDIRk, HDIRq, p, k, q, mu1_p, mu1_k, mu1_q, pj, kj, qj, nu, t);
Sp_raw = 0.5 * (term1_p + term2_p);

% q-leg: Sq_raw = 0.5*(kernel3(q,k,p) + kernel3(q,p,k))
term1_q = kernel3(E0k, E0p, E0q, HPOLk, HPOLp, HPOLq, HDIRk, HDIRp, q, k, p, mu1_q, mu1_k, mu1_p, qj, kj, pj, nu, t);
term2_q = kernel3(E0p, E0k, E0q, HPOLp, HPOLk, HPOLq, HDIRp, HDIRk, q, p, k, mu1_q, mu1_p, mu1_k, qj, pj, kj, nu, t);
Sq_raw = 0.5 * (term1_q + term2_q);

delta = (Sk_raw + Sp_raw + Sq_raw) / 3.0;
dEk_pol = (Sk_raw - delta) * dv;
dEp_pol = (Sp_raw - delta) * dv;
dEq_pol = (Sq_raw - delta) * dv;
s = dEk_pol + dEp_pol + dEq_pol;
dEk_pol = dEk_pol - s/3.0;
dEp_pol = dEp_pol - s/3.0;
dEq_pol = dEq_pol - s/3.0;

% ============================================================================
% 4. ST_NL_ISO: kernel4
% ============================================================================
% k-leg: STk_raw = 0.5*(kernel4(k,p,q) + kernel4(k,q,p))
term1_k = kernel4(Ep, E0Tq, E0Tk, k, p, q, mu3_k, nu, D, t-t0);
term2_k = kernel4(Eq, E0Tp, E0Tk, k, q, p, mu3_k, nu, D, t-t0);
STk_raw = 0.5 * (term1_k + term2_k);

% p-leg: STp_raw = 0.5*(kernel4(p,q,k) + kernel4(p,k,q))
term1_p = kernel4(Eq, E0Tk, E0Tp, p, q, k, mu3_p, nu, D, t-t0);
term2_p = kernel4(Ek, E0Tq, E0Tp, p, k, q, mu3_p, nu, D, t-t0);
STp_raw = 0.5 * (term1_p + term2_p);

% q-leg: STq_raw = 0.5*(kernel4(q,k,p) + kernel4(q,p,k))
term1_q = kernel4(Ek, E0Tp, E0Tq, q, k, p, mu3_q, nu, D, t-t0);
term2_q = kernel4(Ep, E0Tk, E0Tq, q, p, k, mu3_q, nu, D, t-t0);
STq_raw = 0.5 * (term1_q + term2_q);

delta = (STk_raw + STp_raw + STq_raw) / 3.0;
dETk_iso = (STk_raw - delta) * dv;
dETp_iso = (STp_raw - delta) * dv;
dETq_iso = (STq_raw - delta) * dv;
s = dETk_iso + dETp_iso + dETq_iso;
dETk_iso = dETk_iso - s/3.0;
dETp_iso = dETp_iso - s/3.0;
dETq_iso = dETq_iso - s/3.0;

% ============================================================================
% 5. ST_NL_DIR: kernel5
% ============================================================================
% k-leg: STk_raw = 0.5*(kernel5(k,p,q) + kernel5(k,q,p))
term1_k = kernel5(E0p, E0q, E0Tk, E0Tp, E0Tq, HPOLp, HPOLq, HDIRp, HDIRq, HTk, HTp, HTq, k, p, q, mu3_k, nu, D, t-t0);
term2_k = kernel5(E0q, E0p, E0Tk, E0Tq, E0Tp, HPOLq, HPOLp, HDIRq, HDIRp, HTk, HTq, HTp, k, q, p, mu3_k, nu, D, t-t0);
STk_raw = 0.5 * (term1_k + term2_k);

% p-leg: STp_raw = 0.5*(kernel5(p,q,k) + kernel5(p,k,q))
term1_p = kernel5(E0q, E0k, E0Tp, E0Tq, E0Tk, HPOLq, HPOLk, HDIRq, HDIRk, HTp, HTq, HTk, p, q, k, mu3_p, nu, D, t-t0);
term2_p = kernel5(E0k, E0q, E0Tp, E0Tk, E0Tq, HPOLk, HPOLq, HDIRk, HDIRq, HTp, HTk, HTq, p, k, q, mu3_p, nu, D, t-t0);
STp_raw = 0.5 * (term1_p + term2_p);

% q-leg: STq_raw = 0.5*(kernel5(q,k,p) + kernel5(q,p,k))
term1_q = kernel5(E0k, E0p, E0Tq, E0Tk, E0Tp, HPOLk, HPOLp, HDIRk, HDIRp, HTq, HTk, HTp, q, k, p, mu3_q, nu, D, t-t0);
term2_q = kernel5(E0p, E0k, E0Tq, E0Tp, E0Tk, HPOLp, HPOLk, HDIRp, HDIRk, HTq, HTp, HTk, q, p, k, mu3_q, nu, D, t-t0);
STq_raw = 0.5 * (term1_q + term2_q);

delta = (STk_raw + STp_raw + STq_raw) / 3.0;
dETk_dir = (STk_raw - delta) * dv;
dETp_dir = (STp_raw - delta) * dv;
dETq_dir = (STq_raw - delta) * dv;
s = dETk_dir + dETp_dir + dETq_dir;
dETk_dir = dETk_dir - s/3.0;
dETp_dir = dETp_dir - s/3.0;
dETq_dir = dETq_dir - s/3.0;

% ============================================================================
% 6. SF_NL: kernel6
% ============================================================================
% ALL arguments fully permuted for each leg
% k-leg: SFk_raw = 0.5*(kernel6(k,p,q) + kernel6(k,q,p))
term1_k = kernel6(E0p, E0q, E0k, EFp, EFq, EFk, k, p, q, mu3_k, mu3_p, mu3_q, nu, D, t-t0);
term2_k = kernel6(E0q, E0p, E0k, EFq, EFp, EFk, k, q, p, mu3_k, mu3_q, mu3_p, nu, D, t-t0);
SFk_raw = 0.5 * (term1_k + term2_k);

% p-leg: SFp_raw = 0.5*(kernel6(p,q,k) + kernel6(p,k,q))
term1_p = kernel6(E0q, E0k, E0p, EFq, EFk, EFp, p, q, k, mu3_p, mu3_q, mu3_k, nu, D, t-t0);
term2_p = kernel6(E0k, E0q, E0p, EFk, EFq, EFp, p, k, q, mu3_p, mu3_k, mu3_q, nu, D, t-t0);
SFp_raw = 0.5 * (term1_p + term2_p);

% q-leg: SFq_raw = 0.5*(kernel6(q,k,p) + kernel6(q,p,k))
term1_q = kernel6(E0k, E0p, E0q, EFk, EFp, EFq, q, k, p, mu3_q, mu3_k, mu3_p, nu, D, t-t0);
term2_q = kernel6(E0p, E0k, E0q, EFp, EFk, EFq, q, p, k, mu3_q, mu3_p, mu3_k, nu, D, t-t0);
SFq_raw = 0.5 * (term1_q + term2_q);

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
% KERNEL FUNCTIONS
% ============================================================================

% kernel1: S_NL_ISO (Isotropic energy transfer)
function val = kernel1(E0_p, E0_q, E0_k, k, p, q, mu1_k, mu1_p, mu1_q, kj, pj, qj, nu, t)
x = (k^2 + p^2 - q^2) / (2*k*p);
y = (k^2 + q^2 - p^2) / (2*k*q);
z = (p^2 + q^2 - k^2) / (2*p*q);
thetaVal = theta(nu, k, p, q, kj, pj, qj, t, mu1_k, mu1_p, mu1_q);
val = thetaVal * 16.0 * pi^2 * p^2 * k^2 * q * (x*y + z^3) * E0_q * (E0_p - E0_k);
end

% kernel2: S_NL_DIR (Directional anisotropy)
function val = kernel2(E0_p, E0_q, E0_k, HPOL_p, HPOL_q, HDIR_p, HDIR_q, HDIR_k, k, p, q, mu1_k, mu1_p, mu1_q, kj, pj, qj, nu, t)
x = (k^2 + p^2 - q^2) / (2*k*p);
y = (k^2 + q^2 - p^2) / (2*k*q);
z = (p^2 + q^2 - k^2) / (2*p*q);
thetaVal = theta(nu, k, p, q, kj, pj, qj, t, mu1_k, mu1_p, mu1_q);

kernel21 = E0_q * (E0_p - E0_k) * HPOL_q;
kernel22 = E0_q * E0_p * HPOL_p;
kernel23 = E0_q * (E0_p - E0_k) * HDIR_q;
kernel24 = E0_q * E0_p * HDIR_p;
kernel25 = E0_q * E0_k * HDIR_k;

val = thetaVal * 4.0 * pi^2 * p^2 * k^2 * q * (...
    (y^2 - 1.0) * (x*y + z^3) * kernel21 + z * (1.0 - z^2)^2 * kernel22) + ...
    thetaVal * 8.0 * pi^2 * p^2 * k^2 * q * (x*y + z^3) * (...
    (3.0*y^2 - 1.0) * kernel23 + (3.0*z^2 - 1.0) * kernel24 - 2.0 * kernel25);
end

% kernel3: S_NL_POL (Poloidal anisotropy)
function val = kernel3(E0_p, E0_q, E0_k, HPOL_p, HPOL_q, HPOL_k, HDIR_p, HDIR_q, k, p, q, mu1_k, mu1_p, mu1_q, kj, pj, qj, nu, t)
x = (k^2 + p^2 - q^2) / (2*k*p);
y = (k^2 + q^2 - p^2) / (2*k*q);
z = (p^2 + q^2 - k^2) / (2*p*q);
thetaVal = theta(nu, k, p, q, kj, pj, qj, t, mu1_k, mu1_p, mu1_q);

kernel31 = E0_q * E0_p * HPOL_p;
kernel32 = E0_q * E0_k * HPOL_k;
kernel33 = E0_q * (E0_p - E0_k) * HPOL_q;
kernel34 = kernel31;
kernel35 = E0_q * E0_k * HPOL_q;
kernel36 = E0_q * (E0_p - E0_k) * HDIR_q;
kernel37 = E0_q * E0_p * HDIR_p;

val = thetaVal * 4.0 * pi^2 * p^2 * k^2 * q * (...
    (x*y + z^3) * ((1.0 + z^2) * kernel31 - 4.0 * kernel32) + ...
    z * (z^2 - 1.0) * (1.0 + y^2) * kernel33 + ...
    2.0 * z * (z^2 - y^2) * kernel34 + ...
    2.0 * y * x * (z^2 - 1.0) * kernel35) + ...
    thetaVal * 24.0 * pi^2 * p^2 * k^2 * q * z * (z^2 - 1.0) * (...
    (y^2 - 1.0) * kernel36 + (z^2 - 1.0) * kernel37);
end

% kernel4: ST_NL_ISO (Scalar isotropic)
function val = kernel4(E_p, E0T_q, E0T_k, k, p, q, mu3_k, nu, D, t_rel)
y = (k^2 + q^2 - p^2) / (2*k*q);
thetaTVal = thetaT(nu, D, k, p, q, t_rel, mu3_k);
kernel4_inner = E_p * (k^2 * E0T_q - p^2 * E0T_k);
val = thetaTVal * k / (p*q) * (1.0 - y^2) * kernel4_inner;
end

% kernel5: ST_NL_DIR (Scalar directional)
function val = kernel5(E0_p, E0_q, E0T_k, E0T_p, E0T_q, HPOL_p, HPOL_q, HDIR_p, HDIR_q, HT_k, HT_p, HT_q, k, p, q, mu3_k, nu, D, t_rel)
x = (k^2 + p^2 - q^2) / (2*k*p);
y = (k^2 + q^2 - p^2) / (2*k*q);
z = (p^2 + q^2 - k^2) / (2*p*q);
thetaTVal = thetaT(nu, D, k, p, q, t_rel, mu3_k);

kernel51 = E0_q * (E0T_p - E0T_k) * HPOL_q;
kernel52 = E0_q * (E0T_p - E0T_k) * HDIR_q;
kernel53 = E0_q * E0T_p * HT_p;
kernel54 = E0_q * 2.0 * E0T_k * HT_k;

val = 4.0 * thetaTVal * pi^2 * k^2 * p^2 * q * (x*y + z) * (y^2 - 1.0) * kernel51 + ...
    8.0 * thetaTVal * pi^2 * k^2 * p^2 * q * (x*y + z) * (3.0*y^2 - 1.0) * kernel52 + ...
    8.0 * thetaTVal * pi^2 * k^2 * p^2 * q * (x*y + z) * ((3.0*z^2 - 1.0) * kernel53 - kernel54);
end

% kernel6: SF_NL (Flux transfer)
function val = kernel6(E0_p, E0_q, E0_k, EF_p, EF_q, EF_k, k, p, q, mu3_k, mu3_p, mu3_q, nu, D, t_rel)
x = (k^2 + p^2 - q^2) / (2*k*p);
y = (k^2 + q^2 - p^2) / (2*k*q);
z = (p^2 + q^2 - k^2) / (2*p*q);

thetaF_1 = thetaF(nu, D, k, p, q, t_rel, mu3_p, mu3_q);
thetaF_2 = thetaF(nu, D, p, k, q, t_rel, mu3_k, mu3_q);

kernel61 = E0_q * EF_p;
kernel62 = E0_q * EF_k;
kernel63 = E0_k * EF_q;
kernel64 = E0_k * EF_p;
kernel65 = E0_p * EF_q;
kernel66 = E0_p * EF_k;

val = 4.0 * pi^2 * thetaF_1 * k^2 * p * q * (...
    k * kernel61 * (1.0 + y^2 - z^2 - x*y*z - 2.0*y^2*z^2) - ...
    2.0 * q * (y^3 + x*z) * kernel62) + ...
    4.0 * pi^2 * thetaF_2 * k^2 * p * q * (...
    q * z * (2.0*x*y^2 + y*z - x) * kernel63 - ...
    p * y * (x + y*z) * kernel64 + ...
    k * ((1.0 - y^2 + z^2 - x*y*z - 2.0*y^2*z^2) * kernel65 - 2.0*(1.0 - y^2) * kernel66));
end

% ============================================================================
% HELPER FUNCTIONS
% ============================================================================

function val = interp_log(k_eval, logk, logval)
logk_eval = log(k_eval);
logval_interp = interp1(logk, logval, logk_eval, 'linear', 'extrap');
val = exp(logval_interp);
end

function val = interp_log_signed(k_eval, logk, logval, signval)
logk_eval = log(k_eval);
sign_interp = interp1(logk, signval, logk_eval, 'linear', 'extrap');
logval_interp = interp1(logk, logval, logk_eval, 'linear', 'extrap');
val = sign_interp * exp(logval_interp);
end

function [sum_new, c_new] = kahan_add(sum_old, c_old, x)
y = x - c_old;
t = sum_old + y;
c_new = (t - sum_old) - y;
sum_new = t;
end
