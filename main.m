L = 10.0; %lengthscale of the largest eddies
Re_bTARGET = 100.0;
deltaTstarET_Target = 200.0; %how long to run ET after ET is initialized.





%begin local code



termFlag = 0;
NotintegerVal = 0.0;
    
     

     Re_l = 1.0e7;
  

     eta = Re_l^(-3.0/4.0)* L;
     nu = 0.05;
     
     epzilon = nu^3.0/(eta^4.0);
    
     t = 0.0;
     sigma = 2.0;
    
     kEta = 1./eta;
    
     
   
     Pr = 200.0;
     D = nu/Pr;

     
     K0 = 0.01; %manual choice will give 148 kVals
     du = 10.0^(1.0/17.0);

     kVals(1) = K0;
     K0minus1 = K0/du;
     
     kValsExtended(1)=K0/du;
     kValsExtended(2) = K0;
     y = 2; 
       while (kVals(y-1)< (10.0*kEta*sqrt(Pr) ) )
        %if(kVals(y-1)<100)||(kVals(y-1)>kEta*sqrt(Pr))
        kVals(y) = du*kVals(y-1);
        %else
       %     kVals(y) = kVals(y-1)+(kVals(y-1)-kVals(y-2));
      %  end
        kValsExtended(y+1) = kVals(y);
        y = y+1;
       end
       kValsExtended(y+1) = du*kVals(y-1); 

     lengthKvals = length(kVals)
     kLength = length(kVals)
    % [qjMin,qjMax] = triadBoundaries(kVals,lengthKvals); %finds the qj values for kj and pj that satisfy domain condition
 for kj = 1:kLength
        k = kVals(kj);

        fl(kj) = ( k*L/ (((k*L)^1.5+1.5-sigma/4.0)^(2.0/3.0)) )^(5.0/3.0+sigma);
        fEta(kj) = exp(-5.3*(((k*eta)^4.0+0.4^4.0)^(1.0/4.0)-0.4));

        %E= K0*k.^(-5/3)*epzilon^(2/3).*fl(k).*fEta(k*eta);
        E(kj)=  k^(-5.0/3.0) *fl(kj)*fEta(kj);

 end 

     c01 = epzilon^(1./3.) / (2.*nu*trapezoidalIntegration(kVals,kVals.^2.*E));
     c0 = epzilon^(1./3.) / (2.*nu*trapz(kVals,kVals.^2.*E));
     epsTest = 2*nu*trapezoidalIntegration(kVals,kVals.^2.*E*epzilon^(2./3.)*c0);

      E = c0*epzilon^(2./3.)*E;

  KE0 = 0.5*trapezoidalIntegration(kVals,E);

     TauL = KE0/epzilon; %eddy turnover time

     TauEta = (nu/epzilon)^0.5 ;
     TauL2 = TauEta*Re_l^0.5;
     TauL0 = TauL;
 NtauInv = 0.0 ; 
 Re0 = KE0^2/(epzilon*nu)


Eic = icGenerate(kVals);
loglog(kVals,Eic);
%mu = getMu(Eic,kVals);
%mu1 = A1*mu;
%[weight,ySquaredAvg,triadFlag,outsideCutCell,insideCutCell,CxVals,CyVals,Q11]=midpoint2dShoelace(kVals);
kmin = 0.01;
[weight, CxVals,CyVals, dv, edges] = buildTriadWeightsCentroidsExact(kVals, kmin);

 symErr = max(abs(weight - permute(weight,[2 1 3])),[],'all');
fprintf('symErr=%g\n', symErr);

A_k = zeros(length(kVals),1);
for kj=1:length(kVals)
    A_k(kj) = sum(weight(:,:,kj),'all');
end
fprintf('A_k min/max = %g  %g\n', min(A_k), max(A_k));


[S_NL_E, diag] = transfer_scatter_add_kernelE0(kVals, edges, E, weight, CxVals, CyVals);
fprintf('FV total energy transfer = %.20e\n', diag.FV_total_energy_transfer);
fprintf('max triad energy residual = %.20e\n', diag.maxTriadEnergyResidual);
fprintf('triads used = %d\n', diag.numTriadsUsed);

partial = cumsum(S_NL_E .* diff(edges));
semilogx(sqrt(edges(1:end-1).*edges(2:end)), partial), grid on, yline(0,'--')
title('FV partial sums of energy transfer (scatter-add, your kernel)')


dk   = diff(edges(:));
