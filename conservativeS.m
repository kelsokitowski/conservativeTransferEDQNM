function [S_NL_ISO, maxResidual, Tkj] = conservativeS(kVals, edges, kj, E, CxVals, CyVals, weight)
% conservativeS - Compute nonlinear transfer using midpoint rule
%
% Uses area-based weights from midpoint2dShoelace and centroids to evaluate
% the EDQNM triad kernel at cell centers, implementing the classic midpoint rule.
%
% INPUTS:
%   kVals     - wavenumber grid centers [N×1]
%   edges     - wavenumber bin edges [N+1×1]
%   kj        - index for current k wavenumber
%   E         - energy spectrum E(k) [N×1]
%   CxVals    - p-centroids [N×N×N] from centroidOfEdgeCell
%   CyVals    - q-centroids [N×N×N] from centroidOfEdgeCell
%   weight    - area weights [N×N×N] from midpoint2dShoelace
%
% OUTPUTS:
%   S_NL_ISO  - nonlinear transfer rate at k = kVals(kj)
%   maxResidual - maximum conservation residual for this k
%   Tkj       - triad contribution array [N×N×N] for debugging

    N = length(kVals);
    k = kVals(kj);
    dk = edges(kj+1) - edges(kj);

    % Initialize outputs
    S_NL_ISO = 0.0;
    Tkj = zeros(N, N, N);
    maxResidual = 0.0;

    % Spectral density conversion: E0 = E/(4π k²)
    E0 = E ./ (4*pi * kVals.^2);

    % Loop over all (p, q) pairs with non-zero weight
    for pj = 1:N
        for qj = 1:N
            % Skip if no weight (no contribution to integral)
            if weight(pj, qj, kj) == 0
                continue;
            end

            % Get centroids where kernel will be evaluated
            p = CxVals(pj, qj, kj);
            q = CyVals(pj, qj, kj);

            % Validate centroids
            if p == 0 || q == 0
                continue; % Skip invalid centroids
            end

            % Check triad condition at centroid
            if ~triadCondition(q, p, k)
                fprintf('WARNING: centroid (p=%.6e, q=%.6e) violates triad for k=%.6e at (pj=%d,qj=%d,kj=%d)\n', ...
                        p, q, k, pj, qj, kj);
                continue;
            end

            % Compute direction cosines at centroid (p, q, k)
            [x, y, z] = waveCosines(k, p, q);

            % EDQNM kernel evaluation at centroid
            % This is the core physics: transfer rate density in (p,q) space

            % Spectral densities at bin centers
            E0_k = E0(kj);
            E0_p = E0(pj);
            E0_q = E0(qj);

            % Geometric factor (simplified - assumes inviscid, no time decorrelation)
            % Full version would include theta(k,p,q,t) time decorrelation
            b_kpq = (p*q/(k^2)) * y * z;

            % EDQNM transfer kernel (simplified isotropic version)
            % S_k = ∫∫ b(k,p,q) * [p E0_p E0_q - k E0_k (E0_p + E0_q)] dp dq
            kernel_k = b_kpq * (p * E0_p * E0_q - k * E0_k * (E0_p + E0_q));

            % Area weight for midpoint rule: ∫∫ f(p,q) dp dq ≈ Σ f(p_c, q_c) * area
            % Note: weight already includes the area, so no additional dk_p or dk_q
            dS_k = kernel_k * weight(pj, qj, kj);

            % Accumulate contribution
            S_NL_ISO = S_NL_ISO + dS_k;

            % Store for energy conservation check
            Tkj(kj, pj, qj) = dS_k / dk;  % Normalize by dk for triad check

            % ENERGY CONSERVATION CHECK
            % For each triad, we should have: T_k + T_p + T_q = 0
            % We need to check cyclic permutations if all three are computed

            % Find indices where we can compute the partner triads
            % (This requires reverse lookup from centroid to index - simplified here)

            % Simplified residual check: just accumulate magnitude
            % Full check would require storing all three cyclic permutations
            residual = abs(dS_k); % Placeholder - true check needs all 3 terms
            maxResidual = max(maxResidual, residual);
        end
    end

    % Report statistics
    fprintf('k(kj=%d) = %.6e: S_NL_ISO = %.6e, maxResidual = %.6e\n', ...
            kj, k, S_NL_ISO, maxResidual);
end
