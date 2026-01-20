  function [x,y,z]= waveCosines(k,p,q)
      
      %cosines of interior angles of triads
      x = -(k^2.-p^2.-q^2.)/(2.*p*q);
      y = -(p^2.-k^2.-q^2.)/(2.*k*q);
      z = -(q^2.-k^2.-p^2.)/(2.*k*p);
    end 