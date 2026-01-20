
function [A,A2,Cx,Cy] = centroidEdge2(xVals,yVals)
% Assume xVals, yVals define a closed polygon (last point = first point)

% --- 1. Center & scale coordinates for numerical stability ---
xShift = mean(xVals);
yShift = mean(yVals);

xScaled = xVals - xShift;
yScaled = yVals - yShift;

scale = max(abs([xScaled(:); yScaled(:)]));
if scale == 0
    A = 0; Cx = xVals(1); Cy = yVals(1);
    return
end
xScaled = xScaled / scale;
yScaled = yScaled / scale;

% --- 2. Shoelace area with Kahan compensation ---
A = 0.0; cA = 0.0;
for i = 1:(length(xScaled)-1)
    term = xScaled(i)*yScaled(i+1) - xScaled(i+1)*yScaled(i);
    y = term - cA;
    t = A + y;
    cA = (t - A) - y;
    A = t;
end
A = 0.5 * A;      % area in scaled units

% --- 3. Centroid components (also in scaled units) ---
Cx = 0.0; Cy = 0.0;
for i = 1:(length(xScaled)-1)
    cross = xScaled(i)*yScaled(i+1) - xScaled(i+1)*yScaled(i);
    Cx = Cx + (xScaled(i) + xScaled(i+1)) * cross;
    Cy = Cy + (yScaled(i) + yScaled(i+1)) * cross;
end
Cx = Cx / (6.0 * A);
Cy = Cy / (6.0 * A);

% --- 4. Rescale to original coordinates ---
A  = A * scale^2;          % restore area units
Cx = Cx * scale + xShift;  % move centroid back to original location
Cy = Cy * scale + yShift;

% --- 5. (Optional) alternate area check using trapezoid form ---
A2 = 0.0; cA2 = 0.0;
for i = 1:(length(xScaled)-1)
    term = (yScaled(i)+yScaled(i+1))*(xScaled(i)-xScaled(i+1));
    y = term - cA2;
    t = A2 + y;
    cA2 = (t - A2) - y;
    A2 = t;
end
A2 = 0.5 * A2 * scale^2;

relTol = 1e-8;
if abs(A2 - A) > 100.0*relTol * min(abs(A), abs(A2))
    warning('uhoh centroidOfEdgeCell: A2 ≠ A (difference exceeds tolerance)')
end

% --- 6. Results ---
%fprintf('Area = %.12g\n', A);
%fprintf('Centroid = (%.12g, %.12g)\n', Cx, Cy);
end
