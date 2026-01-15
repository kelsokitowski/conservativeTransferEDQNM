 function thetaVal = theta(nu,k,p,q,kj,pj,qj,t,mu1_k,mu1_p,mu1_q)
      
      thetaVal = (1.-exp(-(nu*(k^2.+p^2.+q^2.) + mu1_k+mu1_p+mu1_q )*t) )/ ...
     (nu*(k^2.+p^2.+q^2.) + mu1_k+mu1_p+mu1_q);





     %   thetaVal = 1.0/ ...
     % (nu*(k^2.+p^2.+q^2.) + mu1(kj)+mu1(pj)+mu1(qj) );

     %thetaVal = 1; %debug test

   end
