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


%Eic = icGenerate(kVals);
%loglog(kVals,Eic);
%mu = getMu(Eic,kVals);
%mu1 = A1*mu;
%[weight,ySquaredAvg,triadFlag,outsideCutCell,insideCutCell,CxVals,CyVals,Q11]=midpoint2dShoelace(kVals);
kmin = 0.01;
[weight_k, CxVals_k, CyVals_k, weight_p, CxVals_p, CyVals_p, weight_q, CxVals_q, CyVals_q, dv, edges] = buildTriadWeightsCentroidsExact(kVals, kmin);

% Check cyclic symmetry of dv: dv(i,j,k) = dv(j,k,i) = dv(k,i,j)
cycErr1 = max(abs(dv - permute(dv,[2 3 1])),[],'all');  % (i,j,k) vs (j,k,i)
cycErr2 = max(abs(dv - permute(dv,[3 1 2])),[],'all');  % (i,j,k) vs (k,i,j)
fprintf('Cyclic symmetry errors in dv: (i,j,k)→(j,k,i): %g  (i,j,k)→(k,i,j): %g\n', cycErr1, cycErr2);

% Check p<->q symmetry: dv(i,j,k) = dv(j,i,k)
pqErr = max(abs(dv - permute(dv,[2 1 3])),[],'all');
fprintf('p<->q symmetry error in dv: %g\n', pqErr);

% Check p<->q symmetry for each weight set
symErr_k = max(abs(weight_k - permute(weight_k,[2 1 3])),[],'all');
symErr_p = max(abs(weight_p - permute(weight_p,[2 1 3])),[],'all');
symErr_q = max(abs(weight_q - permute(weight_q,[2 1 3])),[],'all');
fprintf('p<->q symmetry errors in weights: k-slice=%g  p-slice=%g  q-slice=%g\n', symErr_k, symErr_p, symErr_q);

% Verify cyclic consistency of weights: all should equal dv/dk(3rd index)
% Since all three use same formula, they should be identical!
weightErr_kp = max(abs(weight_k - weight_p),[],'all');
weightErr_kq = max(abs(weight_k - weight_q),[],'all');
fprintf('Weight consistency: max|weight_k - weight_p|=%g  max|weight_k - weight_q|=%g\n', weightErr_kp, weightErr_kq);

% Check total area for k-slices (should be consistent across all k)
A_k = zeros(length(kVals),1);
for kj=1:length(kVals)
    A_k(kj) = sum(weight_k(:,:,kj),'all');
end
fprintf('A_k min/max = %g  %g\n', min(A_k), max(A_k));

[S_NL_E, FV_total_energy_transfer, maxTriadEnergyResidual, numTriadsUsed] = transfer_scatter_add_kernelE0(kVals, edges, E, weight_k, CxVals_k, CyVals_k, weight_p, CxVals_p, CyVals_p, weight_q, CxVals_q, CyVals_q);
fprintf('FV total energy transfer = %.20e\n', FV_total_energy_transfer);
fprintf('max triad energy residual = %.20e\n', maxTriadEnergyResidual);
fprintf('triads used = %d\n', numTriadsUsed);

partial = cumsum(S_NL_E .* diff(edges));
semilogx(sqrt(edges(1:end-1).*edges(2:end)), partial), grid on, yline(0,'--')
title('FV partial sums of energy transfer (scatter-add, your kernel)')


dk   = diff(edges(:));
