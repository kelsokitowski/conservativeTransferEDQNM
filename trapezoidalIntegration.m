function theResult = trapezoidalIntegration(xvals,Fvals)
    
     theResult = 0.;
xiRange = length(xvals);

     for j = 2:length(xvals)
        h = (xvals(j)-xvals(j-1))/2.;
        theResult = theResult + h*(Fvals(j)+Fvals(j-1));
     end 
    % trapezoidalIntegration = theResult;

   end 