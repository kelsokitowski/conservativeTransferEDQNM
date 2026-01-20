function [pM,pP,qM,qP] =getCorners(pj,qj,kLength,p,q,kVals)
% getCorners - TEMPORARILY REVERTED to arithmetic mean for debugging
%
% This is the ORIGINAL version using arithmetic mean spacing.
% Need to verify this gives good accuracy before switching to geometric.

if (pj<kLength)

    pP = p+(kVals(pj+1)-kVals(pj))/2;
else
    pP = p;
end
if (pj>1)
    pM = p-(p-kVals(pj-1))/2;
else
    pM = p;

end
if (qj<kLength)
    qP = q+(kVals(qj+1)-q)/2;
else
    qP = q;
end
if (qj>1)
    qM = q-(q-kVals(qj-1))/2;
else
    qM = q;

end


end
