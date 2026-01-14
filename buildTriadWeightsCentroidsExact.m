function [weight, centroidX, centroidY, dv, edges] = buildTriadWeightsCentroidsExact(kVals, kmin)
% buildTriadWeightsCentroidsExact (GEOMETRIC centroid version)
%
% Same dv and weights as before.
% The only change is that the returned evaluation point (centroidX/Y) is now
% the GEOMETRIC centroid (pc,qc) of the dominant convex polygon piece, which
% is guaranteed to lie inside that piece (and hence inside the triad domain
% intersection in the cell). This avoids the failure mode of L-weighted
% centroids (mp_loc/dv_loc) drifting outside in thin slivers.

    if nargin < 2 || isempty(kmin)
        kmin = 0.01;
    end

    kVals = kVals(:);
    N = length(kVals);

    if any(kVals <= 0) || any(~isfinite(kVals))
        error('kVals must be positive finite centers.');
    end

    % -------------------------------
    % 1) Construct log-bin edges
    % -------------------------------
    edges = logCellEdgesFromCenters_withKmin(kVals, kmin);
    dk = edges(2:end) - edges(1:end-1);

    % Allocate outputs in (pj,qj,kj)
    dv     = zeros(N,N,N);
    weight = zeros(N,N,N);

    % Choice A (robust for production): initialize centroids to 0
    centroidX = zeros(N,N,N);
    centroidY = zeros(N,N,N);

    % Choice B (debugging): initialize to NaN to catch accidental use
    % centroidX = nan(N,N,N);
    % centroidY = nan(N,N,N);

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
                % 2b) Exact dv and a guaranteed-inside evaluation point
                % ----------------------------------------------------------
                [vol, ~, ~, px_in, qy_in] = dv_moments_oneCell_exact(kL,kU,pL,pU,qL,qU);

                if vol > 0
                    dv(pj,qj,kj)     = vol;
                    weight(pj,qj,kj) = vol / dkj;

                    centroidX(pj,qj,kj) = px_in;
                    centroidY(pj,qj,kj) = qy_in;

                    % Optional: hard sanity checks (enable if debugging)
                    % if ~(isfinite(px_in) && isreal(px_in) && px_in >= pL && px_in <= pU && px_in > 0)
                    %     error('Bad centroidX at (pj,qj,kj)=(%d,%d,%d): %g not in [%g,%g]', pj,qj,kj,px_in,pL,pU);
                    % end
                    % if ~(isfinite(qy_in) && isreal(qy_in) && qy_in >= qL && qy_in <= qU && qy_in > 0)
                    %     error('Bad centroidY at (pj,qj,kj)=(%d,%d,%d): %g not in [%g,%g]', pj,qj,kj,qy_in,qL,qU);
                    % end

                    % Mirror p<->q symmetry
                    if qj ~= pj
                        dv(qj,pj,kj)     = vol;
                        weight(qj,pj,kj) = vol / dkj;

                        centroidX(qj,pj,kj) = qy_in;
                        centroidY(qj,pj,kj) = px_in;
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
% EXACT dv AND MOMENTS FOR ONE CELL (CORE GEOMETRY)
% =========================================================================
function [dv, mp, mq, px_in, qy_in] = dv_moments_oneCell_exact(kL,kU,pL,pU,qL,qU)

    dv = 0.0; mp = 0.0; mq = 0.0;
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

        % classify branches at (pc,qc)
        s    = pc + qc;
        d    = pc - qc;
        absd = abs(d);

        upper = min(kU, s);
        lower = max(kL, absd);

        if upper <= lower
            continue;
        end

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

        if dv_loc <= 0
            continue;
        end

        % total exact integrals
        dv = dv + dv_loc;
        mp = mp + mp_loc;
        mq = mq + mq_loc;

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
