function output = thetaF(nu,D,k,p,q,t,mu3p,mu3q)

output = ( 1.-exp( -(D*k^2.+nu*(p^2. +q^2.) + mu3+ mu3q )*t ) )/...
    ( D*k^2.+nu*(p^2.+q^2.) +mu3q+ mu3p);

end
