function [Cx,Cy,A] =centroidOfEdgeCell(qM,qP,pM,pP,k)

x = [];
y = [];


%check top cell boundary for intersections with domain boundaries
xVals = [];
yVals = [];

pSolution = qP-k; %top of omega

if ( (pM<pSolution) && (pSolution< pP) )
    
    x = [x, pSolution];
    y = [y, qP];
end
pSolution = qP+k; %bottom right of omega
if ( (pM<pSolution) && (pSolution< pP) )
    x = [x, pSolution];
    y = [y, qP];
end
pSolution = k-qP; %bottom left of omega
if ( (pM<pSolution) && (pSolution< pP) )
    x = [x, pSolution];
    y = [y, qP];
end
if triadCondition(qP,pP,k)
    x = [x, pP];
    y = [y, qP];
end

%sort clockwise
%here y is fixed. Sort x
if ~isempty(x)
    if length(x)>1
xSorted = sort(x,'ascend');
    else
        xSorted = x;
    end
    xVals = [xVals, xSorted];
    yVals = [yVals, y]; 
end


       

x = [];
y = [];

%now do right cell edge (pP fixed)
qSolution = pP-k;
if ( (qM<qSolution) && (qSolution<qP) )
    x = [x,pP];
    y = [y,qSolution];
end
qSolution = pP+k;
if ( (qM<qSolution) && (qSolution<qP) )
    x = [x, pP];
    y = [y, qSolution];
end
qSolution = k-pP;
if ( (qM<qSolution) && (qSolution<qP) )
    x = [x, pP];
    y = [y, qSolution];
end
if triadCondition(qM,pP,k)
    x = [x,pP];
    y = [y,qM];
end

%sort clockwise
%here x = pP is fixed. Sort y decreasing
if ~isempty(x)
    if length(x)>1
ySorted = sort(y,'descend');
    else
        ySorted = y;
    end
    xVals = [xVals, x];
    yVals = [yVals, ySorted]; 
end

x = [];
y = [];

% Now do bottom cell edge (qM fixed)
pSolution = qM-k;
if ( (pM<pSolution) && (pSolution<pP) )
    x = [x, pSolution];
    y = [y, qM];
end
pSolution = qM+k;
if ( (pM<pSolution) && (pSolution<pP) )
    x = [x, pSolution];
    y = [y, qM];
end
if qM>0 %avoid degenerate duplicate (domain degeneracy)
pSolution = k-qM;
if ( (pM<pSolution) && (pSolution<pP) )
    x = [x, pSolution];
    y = [y, qM];
end
end
if triadCondition(qM,pM,k)
    x = [x,pM];
    y = [y,qM];
end
%sort clockwise
%here y=qM is fixed. Sort x descending
if ~isempty(x)
    if length(x)>1
xSorted = sort(x,'descend');
    else
        xSorted = x;
    end
    xVals = [xVals, xSorted];
    yVals = [yVals, y]; 
end

x = [];
y = [];


%now do left cell edge (pM fixed)
qSolution = pM-k;
if ( (qM<qSolution) && (qSolution<qP) )
    x = [x,pM];
    y = [y,qSolution];
end
qSolution = pM+k;
if ( (qM<qSolution) && (qSolution<qP) )
    x = [x, pM];
    y = [y, qSolution];
end
if pM > 0 %avoid duplicate of domain corner
qSolution = k-pM;
if ( (qM<qSolution) && (qSolution<qP) )
    x = [x, pM];
    y = [y, qSolution];
end
end
if triadCondition(qP,pM,k)
    x = [x,pM];
    y = [y,qP];
end

%now sort points in clockwise fashion
%sort clockwise
%here x = pM is fixed. Sort y ascending
if ~isempty(x)
    if length(x)>1
        ySorted = sort(y,'ascend');
    else
        ySorted = y;
    end
    xVals = [xVals, x];
    yVals = [yVals, ySorted]; 
end
if isempty(xVals)
    disp('debug 161 centroidOfEdgeCell.m xVals is empty')
end
xVals = [xVals, xVals(1)];
yVals = [yVals, yVals(1)];%must loop around in gauss formula
xVals = flip(xVals);
yVals = flip(yVals);
%work in local frame
xVals = xVals - pM;
yVals = yVals-qM;
A = 0.0;
for i = 1:(length(xVals)-1)
    A = A + (xVals(i)*yVals(i+1)-xVals(i+1)*yVals(i));
end
A = A/2.0;
A2 = 0.0;
for i = 1:length(xVals)-1
A2 = A2 + (yVals(i)+yVals(i+1))*(xVals(i)-xVals(i+1));
end
A2 = A2/2.0;
if (abs(A2-A)/A*100.0>1e-4)
    disp('uhoh centroidOfEdgeCell A2 neq A ln 181')
end
[A3,A4,Cx2,Cy2] = centroidEdge2(xVals,yVals);
if (abs((A2-A3))>1.0e-05) | (abs((A2-A4))>1.0e-05) | (abs((A-A3))>1.0e-05) | (abs((A-A4))>1.0e-05)
	disp('possible inconsistency with "better" area code')
	fprintf('A = %f , A2 = %f ,A3 = %f , A4 = %f',A,A2,A3,A4)
end
Cx = 0.0;
Cy = 0.0;
for i = 1:(length(xVals)-1)
    Cx = Cx+(xVals(i)+xVals(i+1))*(xVals(i)*yVals(i+1)-xVals(i+1)*yVals(i))/6.0;
    Cy = Cy+(yVals(i)+yVals(i+1))*(xVals(i)*yVals(i+1)-xVals(i+1)*yVals(i))/6.0;
end
%shift back to origin frame.

Cx = Cx/A;
Cy = Cy/A;
A = abs(A); %since clockwise orientation, A will be negative, we want >0
Cx = Cx+pM;
Cy = Cy+qM;
Cx2 = Cx2+pM;
Cy2 = Cy2+qM;
if Cx == 0 || Cy == 0
    disp('centroid at zero?')
end

if abs(Cx-Cx2)> Cx/100.0*1e-04
       disp('centroid x inconsistent result centroidOfEdgeCell :206')
       residual = Cx-Cx2;
       fprintf('Cx-Cx2 = %f, Cx = %f, Cx2 = %f \n',residual,Cx,Cx2) 
end

if abs(Cy-Cy2)> Cy/100.0*1e-04
       disp('centroid y inconsistent result centroidOfEdgeCell :206')
	residual = Cx-Cx2;
       fprintf('Cy-Cy2 = %f, Cy = %f, Cy2 = %f \n',residual,Cy,Cy2)

end
Cy = Cy2; Cx = Cx2;

if (Cy>qP) || (Cx>pP) || (Cx<pM)||(Cy<qM)
    disp('Centroid out of cell. ln 194 centroidOfEdgeCell')
end
if ~triadCondition(Cy,Cx,k)
    disp('centroid does not satisfy triad condition, error')
end


end


