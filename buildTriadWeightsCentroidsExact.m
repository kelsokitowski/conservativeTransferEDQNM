function [weight_k, centroidX_k, centroidY_k, centroidZ_k, weight_p, centroidX_p, centroidY_p, weight_q, centroidX_q, centroidY_q, dv, edges] = buildTriadWeightsCentroidsExact(kVals, kmin)
% buildTriadWeightsCentroidsExact (CYCLIC SYMMETRIC VERSION)
%
% Computes weights and centroids for all three cyclic permutations:
%   (k,p,q), (p,q,k), (q,k,p)
% to achieve energy conservation through geometric symmetry.
%
% Returns:
%   weight_k(pj,qj,kj) = vol/dk(kj) - for k-slice integration (over p,q)
%   weight_p(qj,kj,pj) = vol/dp(pj) - for p-slice integration (over q,k)
%   weight_q(kj,pj,qj) = vol/dq(qj) - for q-slice integration (over k,p)
%   Corresponding centroid arrays for each permutation

    if nargin < 2 || isempty(kmin)
        kmin = 0.01;
    end

    kVals = kVals(:);
    N = length(kVals);

    if any(kVals <= 0) || any(~isfinite(kVals))
        error('kVals must be positive finite centers.');
    end

    % -------------------------------
    % 1) Construct log-bin edges (standard approach)
    % We'll use kVals(kj) directly as k-centroids, not compute from moments
    % -------------------------------
    edges = logCellEdgesFromCenters_withKmin(kVals, kmin);
    dk = edges(2:end) - edges(1:end-1);

    % Allocate outputs for all three cyclic permutations
    dv = zeros(N,N,N);

    % k-slice: integrate over (p,q) for fixed k
    weight_k   = zeros(N,N,N);
    centroidX_k = zeros(N,N,N);  % p-centroid
    centroidY_k = zeros(N,N,N);  % q-centroid
    centroidZ_k = zeros(N,N,N);  % k-centroid

    % p-slice: integrate over (q,k) for fixed p
    weight_p   = zeros(N,N,N);
    centroidX_p = zeros(N,N,N);
    centroidY_p = zeros(N,N,N);

    % q-slice: integrate over (k,p) for fixed q
    weight_q   = zeros(N,N,N);
    centroidX_q = zeros(N,N,N);
    centroidY_q = zeros(N,N,N);

    e = edges;

    % -------------------------------
    % 2) Main loops (k bins outermost)
    % -------------------------------
    for kj = 1:N
        kL = e(kj);
        kU = e(kj+1);
        dkj = dk(kj);

        for pj = 1:N
            pL = e(pj);
            pU = e(pj+1);

            for qj = pj:N
                qL = e(qj);
                qU = e(qj+1);

                % ----------------------------------------------------------
                % 2a) Cheap pruning
                % ----------------------------------------------------------
                if (pU + qU) < kL
                    continue;
                end
                minDiff = intervalMinAbsDiff(pL,pU,qL,qU);
                if minDiff >= kU
                    continue;
                end

                % ----------------------------------------------------------
                % 2b) Exact dv and moments (true centroids)
                % ----------------------------------------------------------
                [vol, mp, mq, mk, ~, ~] = dv_moments_oneCell_exact(kL,kU,pL,pU,qL,qU);

                if vol > 0
                    % Fill dv for all 6 permutations (cyclic + p↔q symmetry)
                    dv(pj,qj,kj) = vol;
                    dv(qj,kj,pj) = vol;  % cyclic: (p,q,k) → (q,k,p)
                    dv(kj,pj,qj) = vol;  % cyclic: (p,q,k) → (k,p,q)
                    if qj ~= pj
                        dv(qj,pj,kj) = vol;  % p↔q mirror
                        dv(pj,kj,qj) = vol;  % p↔q mirror of 2nd cyclic
                        dv(kj,qj,pj) = vol;  % p↔q mirror of 3rd cyclic
                    end

                    % Compute p and q centroids from exact analytical moments (always)
                    pc_true = mp / vol;
                    qc_true = mq / vol;

                    % HYBRID CENTROID STRATEGY FOR k-COORDINATE:
                    % - Interior cells (no boundary cuts): use prescribed kVals(kj)
                    % - Boundary-cut cells: compute from exact analytical moments
                    if is_interior_cell(pL, pU, qL, qU, kL, kU)
                        % Interior cell: use prescribed centroid for consistency
                        kc_true = kVals(kj);
                    else
                        % Boundary-cut cell: compute from exact analytical moments
                        kc_computed = mk / vol;

                        % Validate k-centroid is within physical bounds
                        % Use relative tolerance based on bin width
                        tol = 1e-12 * (kU - kL);
                        if kc_computed >= kL - tol && kc_computed <= kU + tol
                            % Clamp to exact bounds (handles tiny numerical overshoot)
                            kc_true = max(kL, min(kU, kc_computed));
                        else
                            % Truly pathological case: use geometric center fallback
                            % With boundary detection, this should be extremely rare
                            kc_true = 0.5 * (kL + kU);
                            fprintf('WARNING: Invalid k-centroid %.6e outside [%.6e, %.6e] for boundary cell (pj=%d,qj=%d,kj=%d)\n', ...
                                    kc_computed, kL, kU, pj, qj, kj);
                            fprintf('         Using geometric center %.6e instead\n', kc_true);
                        end
                    end

                    % Store centroids for k-slice at (pj,qj,kj)
                    centroidX_k(pj,qj,kj) = pc_true;  % p coordinate
                    centroidY_k(pj,qj,kj) = qc_true;  % q coordinate
                    centroidZ_k(pj,qj,kj) = kc_true;  % k coordinate (hybrid strategy!)

                    % Mirror p<->q symmetry for k-slice centroids
                    if qj ~= pj
                        centroidX_k(qj,pj,kj) = qc_true;
                        centroidY_k(qj,pj,kj) = pc_true;
                        centroidZ_k(qj,pj,kj) = kc_true;  % k-centroid same under p<->q swap
                    end
                end
            end
        end
    end

    % Compute weights from dv (now properly filled with cyclic symmetry)
    % weight_k(pj,qj,kj) = dv(pj,qj,kj) / dk(kj)  -- 3rd index is k
    % weight_p(qj,kj,pj) = dv(qj,kj,pj) / dk(pj)  -- 3rd index is p
    % weight_q(kj,pj,qj) = dv(kj,pj,qj) / dk(qj)  -- 3rd index is q
    for i = 1:N
        for j = 1:N
            for k = 1:N
                if dv(i,j,k) > 0
                    weight_k(i,j,k) = dv(i,j,k) / dk(k);
                    weight_p(i,j,k) = dv(i,j,k) / dk(k);  % All use 3rd index!
                    weight_q(i,j,k) = dv(i,j,k) / dk(k);

                    % Fill bin centers for all locations where weight > 0
                    % (geometric centroids already filled where computed)
                    if centroidX_k(i,j,k) == 0
                        centroidX_k(i,j,k) = kVals(i);
                    end
                    if centroidY_k(i,j,k) == 0
                        centroidY_k(i,j,k) = kVals(j);
                    end
                    if centroidZ_k(i,j,k) == 0
                        centroidZ_k(i,j,k) = kVals(k);
                    end

                    if centroidX_p(i,j,k) == 0
                        centroidX_p(i,j,k) = kVals(i);
                    end
                    if centroidY_p(i,j,k) == 0
                        centroidY_p(i,j,k) = kVals(j);
                    end

                    if centroidX_q(i,j,k) == 0
                        centroidX_q(i,j,k) = kVals(i);
                    end
                    if centroidY_q(i,j,k) == 0
                        centroidY_q(i,j,k) = kVals(j);
                    end
                end
            end
        end
    end
end

% =========================================================================
% EDGE CONSTRUCTION
% =========================================================================
function edges = logCellEdgesFromCenters_withKmin(x, kmin)
    x = x(:);
    N = length(x);

    edges = zeros(N+1,1);
    edges(2:N) = sqrt(x(1:N-1).*x(2:N));

    edges(1) = kmin;
    edges(N+1) = x(N)^2 / edges(N);

    if any(edges <= 0) || any(~isfinite(edges))
        error('Computed edges are nonpositive or nonfinite.');
    end
    if any(diff(edges) <= 0)
        error('Edges must be strictly increasing. Check kmin vs kVals(1).');
    end
end

% =========================================================================
% PRUNING HELPER
% =========================================================================
function md = intervalMinAbsDiff(a,b,c,d)
    if (b >= c) && (d >= a)
        md = 0.0;
    else
        if b < c
            md = c - b;
        else
            md = a - d;
        end
    end
end

% =========================================================================
% BOUNDARY DETECTION FOR HYBRID CENTROID STRATEGY
% =========================================================================
function is_interior = is_interior_cell(pL, pU, qL, qU, kL, kU)
    % Check if cell [pL,pU] × [qL,qU] × [kL,kU] is completely interior
    % to the triad domain |p-q| < k < p+q (no boundary intersections)
    %
    % Interior means: for all (p,q) in [pL,pU]×[qL,qU],
    %                 the entire k-range [kL,kU] is strictly inside the domain
    %
    % Returns true if cell does NOT intersect any boundary surface

    % Use small margin to ensure truly interior (not just touching boundary)
    margin = 1e-10 * max([pU-pL, qU-qL, kU-kL]);

    % Minimum of (p+q) over the box: occurs at (pL, qL)
    min_sum = pL + qL;

    % Maximum of |p-q| over the box: check all corners
    % Since |p-q| increases as we move to opposite corners
    max_diff = max(abs(pU - qL), abs(qU - pL));

    % Interior condition with margin:
    % - kU < min_sum: upper boundary k = p+q is never reached
    % - kL > max_diff: lower boundary k = |p-q| is never reached
    is_interior = (kU < min_sum - margin) && (kL > max_diff + margin);
end

% =========================================================================
% EXACT dv AND MOMENTS FOR ONE CELL (CORE GEOMETRY)
% =========================================================================
function [dv, mp, mq, mk, px_in, qy_in] = dv_moments_oneCell_exact(kL,kU,pL,pU,qL,qU)

    dv = 0.0; mp = 0.0; mq = 0.0; mk = 0.0;
    px_in = 0.0; qy_in = 0.0;

    if kU <= kL
        return;
    end

    rect = [pL qL;
            pU qL;
            pU qU;
            pL qU];

    lines = [
        1  1  kL;
        1  1  kU;
        1 -1  0;    % NEW: split on p = q to handle |p-q| correctly
        1 -1  kL;
        1 -1  kU;
        1 -1 -kL;
        1 -1 -kU
    ];

    polys = {rect};
    for i = 1:size(lines,1)
        a = lines(i,1); b = lines(i,2); c = lines(i,3);

        newPolys = {};
        for t = 1:numel(polys)
            P = polys{t};
            if isempty(P) || size(P,1) < 3
                continue;
            end

            Pin  = clip_halfspace(P,  a,b,c);
            Pout = clip_halfspace(P, -a,-b,-c);

            if ~isempty(Pin)  && size(Pin,1)  >= 3, newPolys{end+1}  = Pin;  end %#ok<AGROW>
            if ~isempty(Pout) && size(Pout,1) >= 3, newPolys{end+1} = Pout; end %#ok<AGROW>
        end

        polys = newPolys;
        if isempty(polys)
            return;
        end
    end

    % Track dominant polygon piece (largest dv_loc) for inside point
    best_dv = 0.0;
    best_px = 0.0;
    best_qy = 0.0;

    for t = 1:numel(polys)
        P = polys{t};
        if isempty(P) || size(P,1) < 3
            continue;
        end

        [A, Mx, My, Ixx, Iyy, Ixy] = poly_moments(P(:,1), P(:,2));
        if A <= 0
            continue;
        end

        % GEOMETRIC centroid of convex polygon piece (guaranteed inside)
        pc = Mx / A;
        qc = My / A;

        % CRITICAL: Check that upper > lower at ALL vertices
        % Since upper and lower are piecewise linear (after splits), they are
        % either constant or linear within each piece, so extrema occur at vertices
        % Use relative tolerance based on bin width
        tol = 1e-10 * (kU - kL);
        valid_piece = true;
        for vi = 1:size(P,1)
            pv = P(vi,1);
            qv = P(vi,2);
            upper_v = min(kU, pv + qv);
            lower_v = max(kL, abs(pv - qv));
            if upper_v <= lower_v + tol
                valid_piece = false;
                break;
            end
        end
        if ~valid_piece
            continue;
        end

        % classify branches at (pc,qc)
        s    = pc + qc;
        d    = pc - qc;
        absd = abs(d);

        upper = min(kU, s);
        lower = max(kL, absd);

        upper_is_kU = (kU <= s + 1e-13);
        lower_is_kL = (kL >= absd - 1e-13);
        d_nonneg    = (d >= 0);

        if upper_is_kU
            au = 0; bu = 0; cu = kU;
        else
            au = 1; bu = 1; cu = 0;
        end

        if lower_is_kL
            al = 0; bl = 0; cl = kL;
        else
            if d_nonneg
                al = 1;  bl = -1; cl = 0;   % p-q
            else
                al = -1; bl =  1; cl = 0;   % q-p
            end
        end

        a = au - al;
        b = bu - bl;
        c = cu - cl;

        % L should be positive at centroid
        Lc = a*pc + b*qc + c;
        if Lc <= 0
            continue;
        end

        dv_loc = a*Mx + b*My + c*A;
        mp_loc = a*Ixx + b*Ixy + c*Mx;
        mq_loc = a*Ixy + b*Iyy + c*My;

        % k-moment: mk = 0.5 * ∫∫ (upper² - lower²) dp dq
        % where upper = au*p + bu*q + cu, lower = al*p + bl*q + cl
        % Using (a²-b²) = (a-b)(a+b) factorization:
        % upper² - lower² = (a*p + b*q + c) * ((au+al)*p + (bu+bl)*q + (cu+cl))
        a_sum = au + al;
        b_sum = bu + bl;
        c_sum = cu + cl;
        mk_loc = 0.5 * (a*a_sum*Ixx + (a*b_sum + b*a_sum)*Ixy + b*b_sum*Iyy ...
                       + (a*c_sum + c*a_sum)*Mx + (b*c_sum + c*b_sum)*My ...
                       + c*c_sum*A);

        if dv_loc <= 0
            continue;
        end

        % total exact integrals
        dv = dv + dv_loc;
        mp = mp + mp_loc;
        mq = mq + mq_loc;
        mk = mk + mk_loc;

        % ---- CRITICAL CHANGE ----
        % Choose geometric centroid as the inside evaluation point.
        % This is guaranteed inside the convex polygon piece and avoids
        % L-weighted centroid drift in sliver/cancellation cases.
        if dv_loc > best_dv
            best_dv = dv_loc;
            best_px = pc;
            best_qy = qc;
        end
    end

    px_in = best_px;
    qy_in = best_qy;

    % Roundoff guard
    if dv < 0 && dv > -1e-12*max(1,abs(dv))
        dv = 0; mp = 0; mq = 0;
        px_in = 0; qy_in = 0;
    end
end

% =========================================================================
% POLYGON CLIPPING: Sutherland–Hodgman
% =========================================================================
function poly2 = clip_halfspace(poly, a, b, c)
    if isempty(poly)
        poly2 = poly; return;
    end

    X = poly(:,1); Y = poly(:,2);
    n = length(X);

    f = a*X + b*Y - c;

    tol = 1e-14;
    poly2 = zeros(0,2);

    for i = 1:n
        i2 = mod(i,n) + 1;

        x1 = X(i);  y1 = Y(i);  f1 = f(i);
        x2 = X(i2); y2 = Y(i2); f2 = f(i2);

        in1 = (f1 <= tol);
        in2 = (f2 <= tol);

        if in1 && in2
            poly2(end+1,:) = [x2 y2]; %#ok<AGROW>
        elseif in1 && ~in2
            t = f1/(f1 - f2);
            poly2(end+1,:) = [x1 + t*(x2-x1), y1 + t*(y2-y1)]; %#ok<AGROW>
        elseif ~in1 && in2
            t = f1/(f1 - f2);
            poly2(end+1,:) = [x1 + t*(x2-x1), y1 + t*(y2-y1)]; %#ok<AGROW>
            poly2(end+1,:) = [x2 y2]; %#ok<AGROW>
        end
    end

    poly2 = cleanup_poly(poly2);
end

function poly = cleanup_poly(poly)
    if isempty(poly), return; end
    tol = 1e-14;

    keep = true(size(poly,1),1);
    for i = 2:size(poly,1)
        if norm(poly(i,:) - poly(i-1,:)) < tol
            keep(i) = false;
        end
    end
    poly = poly(keep,:);

    if size(poly,1) >= 2 && norm(poly(end,:) - poly(1,:)) < tol
        poly(end,:) = [];
    end

    if size(poly,1) < 3
        poly = zeros(0,2);
    end
end

% =========================================================================
% EXACT POLYGON MOMENTS (GREEN'S THEOREM FORMULAS)
% =========================================================================
function [A, Mx, My, Ixx, Iyy, Ixy] = poly_moments(x, y)

    x = x(:); y = y(:);
    n = length(x);

    if n < 3
        A=0; Mx=0; My=0; Ixx=0; Iyy=0; Ixy=0; return;
    end

    x2 = x([2:end 1]);
    y2 = y([2:end 1]);

    cross = x.*y2 - x2.*y;

    A  = 0.5   * sum(cross);
    Mx = (1/6) * sum((x + x2).*cross);
    My = (1/6) * sum((y + y2).*cross);

    Ixx = (1/12) * sum((x.^2 + x.*x2 + x2.^2).*cross);
    Iyy = (1/12) * sum((y.^2 + y.*y2 + y2.^2).*cross);
    Ixy = (1/24) * sum((2*x.*y + x.*y2 + x2.*y + 2*x2.*y2).*cross);

    % Make orientation positive (CCW)
    if A < 0
        A = -A; Mx = -Mx; My = -My;
        Ixx = -Ixx; Iyy = -Iyy; Ixy = -Ixy;
    end
end

