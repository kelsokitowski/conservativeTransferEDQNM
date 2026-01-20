function [pM,pP,qM,qP] =getCorners(pj,qj,kLength,p,q,kVals)
% getCorners - Compute cell edges using GEOMETRIC MEAN spacing
%
% Updated to match buildTriadWeightsCentroidsExact edge definition:
%   edges(i) = sqrt(kVals(i-1) * kVals(i)) for i=2:N
%
% This ensures consistent cell definitions between the two methods.

if (pj < kLength)
    % Upper edge: geometric mean between p and next point
    pP = sqrt(kVals(pj) * kVals(pj+1));
else
    % Last cell: use squared geometric extension
    % Matches: edges(N+1) = kVals(N)^2 / edges(N)
    if pj > 1
        pM_temp = sqrt(kVals(pj-1) * kVals(pj));
        pP = kVals(pj)^2 / pM_temp;
    else
        pP = kVals(pj);  % Single point case
    end
end

if (pj > 1)
    % Lower edge: geometric mean between previous point and p
    pM = sqrt(kVals(pj-1) * kVals(pj));
else
    % First cell: would use kmin, but for consistency use p itself
    pM = kVals(1);  % Could also use kmin parameter if passed in
end

if (qj < kLength)
    % Upper edge: geometric mean between q and next point
    qP = sqrt(kVals(qj) * kVals(qj+1));
else
    % Last cell: use squared geometric extension
    if qj > 1
        qM_temp = sqrt(kVals(qj-1) * kVals(qj));
        qP = kVals(qj)^2 / qM_temp;
    else
        qP = kVals(qj);  % Single point case
    end
end

if (qj > 1)
    % Lower edge: geometric mean between previous point and q
    qM = sqrt(kVals(qj-1) * kVals(qj));
else
    % First cell: would use kmin, but for consistency use q itself
    qM = kVals(1);  % Could also use kmin parameter if passed in
end

end
