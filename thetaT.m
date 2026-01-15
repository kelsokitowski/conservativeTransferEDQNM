 function thetaTval = thetaT(nu,D,k,p,q,t,mu3q)
      

      thetaTval = ( 1.-exp( -(D*(k^2.+p^2.) +nu*q^2. + mu3q )*t ))/ ...
    (D*(k^2.+p^2.) + nu*q^2. +mu3q );

   end

  

  