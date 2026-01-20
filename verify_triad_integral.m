% verify_triad_integral.m
%
% Comprehensive verification of triad integral computation comparing:
% 1. Exact volume integration (buildTriadWeightsCentroidsExact)
% 2. Area-based midpoint rule (midpoint2dShoelace + conservativeS)
%
% Tests energy conservation, spectral transfer accuracy, and midpoint rule validity

clear; close all;

fprintf('=================================================================\n');
fprintf('TRIAD INTEGRAL VERIFICATION - Midpoint Rule vs Exact Integration\n');
fprintf('=================================================================\n\n');

%% Setup test grid
kmin = 0.01;
kmax = 100;
N = 145;  % Grid points

% Generate log-spaced wavenumber grid
kVals = logspace(log10(kmin), log10(kmax), N)';

fprintf('Grid setup:\n');
fprintf('  kmin = %.6e, kmax = %.6e\n', kmin, kmax);
fprintf('  N = %d points\n\n', N);

%% Generate test energy spectrum
fprintf('Generating test energy spectrum...\n');
E = icGenerate(kVals);  % Initial condition spectrum
fprintf('  Total energy = %.6e\n', trapz(kVals, E));
fprintf('  Energy range: [%.6e, %.6e]\n\n', min(E), max(E));

%% Compute weights and centroids using EXACT volume integration
fprintf('Computing exact volume integration (buildTriadWeightsCentroidsExact)...\n');
tic;
[weight_exact, CxVals_exact, CyVals_exact, CzVals_exact, ~, ~, ~, dv_exact, edges] = ...
    buildTriadWeightsCentroidsExact(kVals, kmin);
t_exact = toc;
fprintf('  Time: %.2f seconds\n', t_exact);
fprintf('  Non-zero weights: %d / %d^3 = %.2f%%\n', ...
        nnz(weight_exact), N, 100*nnz(weight_exact)/N^3);

%% Compute weights and centroids using AREA-BASED midpoint rule
fprintf('\nComputing area-based midpoint rule (midpoint2dShoelace)...\n');
tic;
[weight_area, ~, triadFlag, outsideCutCell, insideCutCell, CxVals_area, CyVals_area, Q11] = ...
    midpoint2dShoelace(kVals);
t_area = toc;
fprintf('  Time: %.2f seconds\n', t_area);
fprintf('  Non-zero weights: %d / %d^3 = %.2f%%\n', ...
        nnz(weight_area), N, 100*nnz(weight_area)/N^3);

%% Compare weight arrays
fprintf('\n=================================================================\n');
fprintf('WEIGHT COMPARISON\n');
fprintf('=================================================================\n');

% Element-wise comparison
weight_diff = abs(weight_exact - weight_area);
rel_diff = weight_diff ./ (abs(weight_exact) + 1e-30);

nonzero_exact = (weight_exact > 0);
nonzero_area = (weight_area > 0);
common_nonzero = nonzero_exact & nonzero_area;

fprintf('Weight statistics:\n');
fprintf('  Exact method: %d non-zero weights\n', nnz(nonzero_exact));
fprintf('  Area method:  %d non-zero weights\n', nnz(nonzero_area));
fprintf('  Common:       %d non-zero weights\n', nnz(common_nonzero));
fprintf('  Max absolute difference: %.6e\n', max(weight_diff(:)));
fprintf('  Max relative difference: %.6e\n', max(rel_diff(common_nonzero)));
fprintf('  Mean relative difference: %.6e\n', mean(rel_diff(common_nonzero)));

%% Test simple integral: ∫∫ 1 dA for each k-slice
fprintf('\n=================================================================\n');
fprintf('INTEGRAL TEST: ∫∫ 1 dA (Total area for each k-slice)\n');
fprintf('=================================================================\n');

analytical = @(k) 2*k*max(kVals) - 1.5*k^2;  % Analytical result

res_exact = zeros(N, 1);
res_area = zeros(N, 1);

for kj = 1:N
    k = kVals(kj);
    dk = edges(kj+1) - edges(kj);

    % Exact integration: weight * dk
    integral_exact = sum(sum(weight_exact(:, :, kj))) * dk;
    res_exact(kj) = integral_exact - analytical(k);

    % Area-based midpoint rule: weight already includes area
    integral_area = sum(sum(weight_area(:, :, kj)));
    res_area(kj) = integral_area - analytical(k);
end

fprintf('Exact integration:\n');
fprintf('  Max residual: %.6e\n', max(abs(res_exact)));
fprintf('  Max percent error: %.6e%%\n', 100*max(abs(res_exact./analytical(kVals))));

fprintf('Area-based midpoint:\n');
fprintf('  Max residual: %.6e\n', max(abs(res_area)));
fprintf('  Max percent error: %.6e%%\n', 100*max(abs(res_area./analytical(kVals))));

%% Compute nonlinear transfer using both methods
fprintf('\n=================================================================\n');
fprintf('NONLINEAR TRANSFER COMPUTATION (Simple Test)\n');
fprintf('=================================================================\n');

% Using EXACT integration
fprintf('\nComputing S_NL using exact volume integration...\n');
S_NL_exact = zeros(N, 1);
tic;
for kj = 1:N
    dk = edges(kj+1) - edges(kj);
    k = kVals(kj);

    % Sum over all (p,q) with non-zero weight
    S_k_local = 0;
    for pj = 1:N
        for qj = 1:N
            if weight_exact(pj, qj, kj) > 0
                p = CxVals_exact(pj, qj, kj);
                q = CyVals_exact(pj, qj, kj);

                % Simple test kernel: just use geometric factor
                [x, y, z] = waveCosines(k, p, q);
                kernel = y * z;  % Simplified test

                S_k_local = S_k_local + kernel * weight_exact(pj, qj, kj) * dk;
            end
        end
    end
    S_NL_exact(kj) = S_k_local;
end
t_transfer_exact = toc;
fprintf('  Time: %.2f seconds\n', t_transfer_exact);

% Using AREA-BASED midpoint rule
fprintf('\nComputing S_NL using area-based midpoint rule...\n');
S_NL_area = zeros(N, 1);
tic;
for kj = 1:N
    k = kVals(kj);

    % Sum over all (p,q) with non-zero weight
    S_k_local = 0;
    for pj = 1:N
        for qj = 1:N
            if weight_area(pj, qj, kj) > 0
                p = CxVals_area(pj, qj, kj);
                q = CyVals_area(pj, qj, kj);

                if p == 0 || q == 0
                    continue;
                end

                % Same test kernel
                [x, y, z] = waveCosines(k, p, q);
                kernel = y * z;

                S_k_local = S_k_local + kernel * weight_area(pj, qj, kj);
            end
        end
    end
    S_NL_area(kj) = S_k_local;
end
t_transfer_area = toc;
fprintf('  Time: %.2f seconds\n', t_transfer_area);

%% Compare transfer results
fprintf('\n=================================================================\n');
fprintf('TRANSFER COMPARISON: S_NL(k)\n');
fprintf('=================================================================\n');

S_diff = abs(S_NL_exact - S_NL_area);
S_rel = S_diff ./ (abs(S_NL_exact) + 1e-30);

fprintf('Transfer statistics:\n');
fprintf('  Max |S_exact - S_area|: %.6e\n', max(S_diff));
fprintf('  Max relative error: %.6e%%\n', 100*max(S_rel));
fprintf('  Mean relative error: %.6e%%\n', 100*mean(S_rel));

% Global energy conservation
G_exact = trapz(kVals, S_NL_exact);
G_area = trapz(kVals, S_NL_area);

fprintf('\nGlobal energy conservation:\n');
fprintf('  ∫ S_NL(k) dk (exact):  %.6e\n', G_exact);
fprintf('  ∫ S_NL(k) dk (area):   %.6e\n', G_area);
fprintf('  Difference:            %.6e\n', abs(G_exact - G_area));

%% Visualization
fprintf('\n=================================================================\n');
fprintf('GENERATING PLOTS\n');
fprintf('=================================================================\n');

% Figure 1: Weight comparison
figure('Position', [100 100 1400 600]);

subplot(1,3,1)
kj_test = round(N/3);
imagesc(log10(abs(weight_exact(:,:,kj_test)) + 1e-30));
colorbar;
title(sprintf('Exact weights (k=%g)', kVals(kj_test)));
xlabel('qj'); ylabel('pj');

subplot(1,3,2)
imagesc(log10(abs(weight_area(:,:,kj_test)) + 1e-30));
colorbar;
title(sprintf('Area weights (k=%g)', kVals(kj_test)));
xlabel('qj'); ylabel('pj');

subplot(1,3,3)
imagesc(log10(abs(weight_exact(:,:,kj_test) - weight_area(:,:,kj_test)) + 1e-30));
colorbar;
title('log10|Difference|');
xlabel('qj'); ylabel('pj');

% Figure 2: Integral test residuals
figure('Position', [100 200 1400 500]);

subplot(1,2,1)
loglog(kVals, abs(res_exact), 'o-', 'DisplayName', 'Exact');
hold on;
loglog(kVals, abs(res_area), 's--', 'DisplayName', 'Area-based');
xlabel('k');
ylabel('|Residual|');
title('Integral Test: |∫∫ 1 dA - analytical|');
legend('Location', 'best');
grid on;

subplot(1,2,2)
loglog(kVals, 100*abs(res_exact./analytical(kVals)), 'o-', 'DisplayName', 'Exact');
hold on;
loglog(kVals, 100*abs(res_area./analytical(kVals)), 's--', 'DisplayName', 'Area-based');
xlabel('k');
ylabel('Percent Error (%)');
title('Integral Test: Percent Error');
legend('Location', 'best');
grid on;

% Figure 3: Transfer comparison
figure('Position', [100 300 1400 500]);

subplot(1,2,1)
loglog(kVals, abs(S_NL_exact), 'o-', 'DisplayName', 'Exact');
hold on;
loglog(kVals, abs(S_NL_area), 's--', 'DisplayName', 'Area-based');
xlabel('k');
ylabel('|S_{NL}(k)|');
title('Nonlinear Transfer');
legend('Location', 'best');
grid on;

subplot(1,2,2)
semilogx(kVals, 100*S_rel, 'o-');
xlabel('k');
ylabel('Relative Error (%)');
title('|S_{exact} - S_{area}| / |S_{exact}|');
grid on;

fprintf('\nVerification complete!\n');
fprintf('=================================================================\n');
