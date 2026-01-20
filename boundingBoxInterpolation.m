function [Q11]=boundingBoxInterpolation(Cx,Cy,kj,pj,qj,p,q,Q11)
%this function obtains the  bounding neighbor node index of bottom left corner for using in
%interpolation.

if ( Cx > p ) %then bounding box is to the right

    Q11(1,pj,qj,kj) = pj;
elseif (Cx<p)
    Q11(1,pj,qj,kj) = pj-1;
elseif (abs(Cx-p)<1e-10)
    Q11(1,pj,qj,kj) = pj;
else
    disp('fallthrough case, boundingBoxInterpolation')
end
if Q11(1,pj,qj,kj) == 0
    disp('uhoh Q11(1,pj,qj,kj) = 0 ln 16 boundingBoxInterpolation.m pj = ')
    pj
    qj
    kj
    disp('terminating program')
    exit

end

if ( Cy > q ) %then bounding box is to the right

    Q11(2,pj,qj,kj) = qj;
    if (qj == 133)
       % disp('uhoh boundingBoxInterpolation ln 23')
    end
elseif (Cy<q)
    Q11(2,pj,qj,kj) = qj-1;
elseif (abs(Cy-q)<1e-10)
     Q11(2,pj,qj,kj) = qj;
      if (qj == 133)
        %disp('uhoh boundingBoxInterpolation ln 30')
    end
else
    disp('fallthrough case, boundingBoxInterpolation')
end
if Q11(2,pj,qj,kj) == 0
    disp('uhoh Q11(2,pj... ) ln 30 boundingBoxInterp')
end














end
