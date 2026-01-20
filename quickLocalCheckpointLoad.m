%lets plot
clear
close all

Pr = 200;
bf = 25000;
read_checkpoint_out(bf,Pr,1,1,0)
function S = read_checkpoint_out(bf_value,Pr,plotE_Flag,plotET_Flag,normalizeET_Flag)
    % Read checkpointOut.bin written by Fortran back into MATLAB
    bf = bf_value;
    folder = sprintf('Pr200/param_%g',bf_value);
    %folder = sprintf('param_%g/param',bf_value);
     load(append(folder,'/ESSIC.mat'))
    %folder = sprintf('newRes');
    fname = fullfile(folder, 'icOut.bin');
if ~isfile(fname)
        error('icOut file not found: %s', fname);
    end
    
    fid = fopen(fname, 'r');
    if fid < 0
        error('Cannot open file: %s', fname);
    end
    
    % Read kLength (int32)
    kLength = fread(fid, 1, 'int32');
    fprintf('Reading checkpoint with kLength = %d\n', kLength);
    
    % Read arrays (each is kLength doubles)
    S.v_1 = fread(fid, kLength, 'double');
    S.v_2 = fread(fid, kLength, 'double');
    S.v_3 = fread(fid, kLength, 'double');
    S.v_4 = fread(fid, kLength, 'double');
    S.v_5 = fread(fid, kLength, 'double');
    S.v_6 = fread(fid, kLength, 'double');
    
    % Read scalars
    S.t = fread(fid, 1, 'double');
    S.tStar = fread(fid, 1, 'double');
    S.counter = fread(fid, 1, 'double');
    
    fclose(fid);
    figure()
    Eic = (S.v_1)';
loglog(kVals,Eic)
axis([-inf inf 10^-20 inf])
title('E IC at scalar initialization')
figure()


    fname = fullfile(folder, 'checkpointOut.bin');
   
    if ~isfile(fname)
        error('Checkpoint file not found: %s', fname);
    end
    
    fid = fopen(fname, 'r');
    if fid < 0
        error('Cannot open file: %s', fname);
    end
    
    % Read kLength (int32)
    kLength = fread(fid, 1, 'int32');
    fprintf('Reading checkpoint with kLength = %d\n', kLength);
    
    % Read arrays (each is kLength doubles)
    S.v_1 = fread(fid, kLength, 'double');
    S.v_2 = fread(fid, kLength, 'double');
    S.v_3 = fread(fid, kLength, 'double');
    S.v_4 = fread(fid, kLength, 'double');
    S.v_5 = fread(fid, kLength, 'double');
    S.v_6 = fread(fid, kLength, 'double');
    
    % Read scalars
    S.t = fread(fid, 1, 'double');
    S.tStar = fread(fid, 1, 'double');
    S.counter = fread(fid, 1, 'double');
    
    fclose(fid);
    
    fprintf('Successfully read checkpoint:\n');
    fprintf('  t = %.16e\n', S.t);
    fprintf('  tStar = %.16e\n', S.tStar);
    fprintf('  counter = %d\n', round(S.counter));

     TauL0=S.t/S.tStar; t30 = 30*S.t/S.tStar; Nt = (S.t-t30)*bf

     Enew = (S.v_1)';
     EHdirNew = (S.v_2)';
     EHpolNew = (S.v_3)';
     ETnew = (S.v_4)';
     ETHnew = (S.v_5)';
     Fnew = (S.v_6)';
     tStar = S.tStar;
     t = tStar*TauL0;


    nu = 0.05;
%Pr = 1;
D = nu/Pr;
epzilon = 2*nu*trapz(kVals,kVals.^2.*Enew);
KE = 0.5*trapz(kVals,Enew);
Re = KE^2/epzilon/nu;
Fr = epzilon/KE/bf;
kB = (epzilon/nu/(D^2))^(1/4);
kEta = (nu^3/epzilon)^(-1/4);

if plotE_Flag == 1

loglog(kVals,Enew);
hold on;
end
ETnormalization = 1;
kNormal = 1;
if normalizeET_Flag == 1
    EmaxIndex = find(Enew==max(Enew));
        ETnormalization = ETnew(EmaxIndex);
        kNormal=kVals(EmaxIndex);
end
if plotET_Flag == 1
    
    
loglog(kVals/kNormal,ETnew./ETnormalization,'--g'); 
%loglog(kVals(maxForcingLoc),ETnew(maxForcingLoc),'Xg');
hold on;
%loglog(kVals,ETnew,'-g')


end
kBloc = find(kVals>kB,1);
kEtaLoc = find(kVals>kEta,1);
hold on;
if plotET_Flag
loglog(kVals(kEtaLoc:end)/kNormal,kVals(kEtaLoc:end).^(-1).*ETnew(kEtaLoc)./ETnormalization*kEta,'--m');
end
if plotE_Flag
hold on; loglog(kVals(1:kEtaLoc),kVals(1:kEtaLoc).^-(5/3)*Enew(10),'--')
end
% if plotET_Flag == 1
% hold on;
% loglog(kB/kNormal,ETnew(kBloc)/ETnormalization,'o')
% end
if plotE_Flag == 1
hold on; 
loglog(kEta,Enew(kEtaLoc),'o')
end
kai = D*trapz(kVals,kVals.^2.*ETnew);
gammaVal = kai/epzilon;

%close all;
%legend('Pope Spectrum','E(t)','ET(t)','peak of forcing','k^{-1} line','kB');
%close all;

axis([-inf 10^5.5 10^(-10) inf])
%plot(kVals,Enew./ETnew)
title(['Energies @ t^* = ',num2str(tStar), ' for Pr = ',num2str(Pr),' bf = ',num2str(bf),' \Gamma = ',num2str(gammaVal),' Fr = ',num2str(Fr)] )
if plotET_Flag == 1 && plotE_Flag == 0
       
        leghandle = findall(gcf, 'tag', 'legend');
    legstr = get(leghandle,'String');
    % ensure legstr is a cell, not a string
    if ischar(legstr)
        legstr = mat2cell(legstr); 
    end
    %legstr(end+1) = {' legend(['bf = ',num2str(bf)])'};
    newPiece = append('bf= ',num2str(bf))
    legstr = [legstr newPiece];
    legend([legstr])
end


    if plotET_Flag == 1 && plotE_Flag == 1
    legend('E(t)','ET(t)','k^{-5/3} line','k^{-1} line','kEta','kB');
    end
if plotET_Flag==0 && plotE_Flag == 1 %then want to look at far field dissipation
axis([-inf 10^5 10^(-30) inf])
legend('E(t)','k^-{5/3} line','kEta')
end 

PE = trapz(kVals,ETnew);
KE = trapz(kVals,Enew);
ratio = KE/PE
ratio2=1/ratio


%============calculate individual terms in PDE=============================
E = Enew;
ET = ETnew;
ETH = ETHnew;
EHdir = EHdirNew;
EHpol = EHpolNew;
F=Fnew;
du = 0;
qjMin = 0;
qjMax = 0;
ySquaredAvg = 0;
load(append(folder,'/weightStuff.mat')) 
load(append(folder,'/ExtForcing.mat'))
max(ExtForcing)
1/max(ExtForcing)
     A1 = 0.355;
     A3 = 1.3;

     kLength = length(kVals);
 
mu= getMu(Eic,kVals); %mu3 of 'Spectral modelling for passive scalar dynamics in hom anis turb'
            mu1 = A1*mu;
            mu3 = A3*mu;
[HDIR,HPOL,HT] = getHvals(kLength,Eic,ET,EHdir,EHpol,ETH);
t0 = 0;
tIC = 50*TauL0;
T = zeros([kLength,kLength,kLength]);
test_triad_symmetry(kVals, Eic, mu1, nu, tIC);
% % ---- extract scalar values ----
% kj = 15; pj = 12; qj = 11;
% k = kVals(kj);
% p = kVals(pj);
% q = kVals(qj);
% 
% fprintf("k=%.15e  p=%.15e  q=%.15e\n", k,p,q);
% 
% % ---- compute direction cosines ----
% [x_kpq , y_kpq , z_kpq] = waveCosines(k,p,q);
% [x_pqk , y_pqk , z_pqk] = waveCosines(p,q,k);
% [x_qkp , y_qkp , z_qkp] = waveCosines(q,k,p);
% 
% fprintf("cosines  (k|p,q):  x=%.15e  y=%.15e  z=%.15e\n", x_kpq,y_kpq,z_kpq);
% fprintf("cosines  (p|q,k):  x=%.15e  y=%.15e  z=%.15e\n", x_pqk,y_pqk,z_pqk);
% fprintf("cosines  (q|k,p):  x=%.15e  y=%.15e  z=%.15e\n", x_qkp,y_qkp,z_qkp);
% 
% % ---- energies & mu ----
% E0_k = E(kj)/(k^2)/(4*pi);
% E0_p = E(pj)/(p^2)/(4*pi);
% E0_q = E(qj)/(q^2)/(4*pi);
% 
% mu_k = mu1(kj);
% mu_p = mu1(pj);
% mu_q = mu1(qj);
% 
% fprintf("E0(k)=%.15e  E0(p)=%.15e  E0(q)=%.15e\n", E0_k,E0_p,E0_q);
% fprintf("mu(k)=%.15e  mu(p)=%.15e  mu(q)=%.15e\n", mu_k,mu_p,mu_q);
% 
% % ---- kernel contributions ----
% Ak = dT_leg(kj,pj,qj,kVals,Eic,mu1,nu,t);
% Ap = dT_leg(pj,qj,kj,kVals,Eic,mu1,nu,t);
% Aq = dT_leg(qj,kj,pj,kVals,Eic,mu1,nu,t);
% 
% fprintf("Ak=%.15e  Ap=%.15e  Aq=%.15e\n", Ak,Ap,Aq);
% fprintf("Ak+Ap+Aq = %.15e\n", Ak+Ap+Aq);
% S_NL_ISO_debug = zeros(size(kVals));
% E_test = kVals;
% for kj = 1:kLength
%     [S_iso,S_NL_DIR,S_NL_POL,ST_NL_ISO,ST_NL_DIR,SF_NL,phi33,Tkj,S_interior,S_cut,S_int] = getTheIndividualTerms(0,nu,D,kj,kVals(kj),kVals,HPOL,HDIR,HT,...
%                                                   E_test,F,ET,mu1,mu3,tIC,weight,...
%                                                   triadFlag,outsideCutCell,Q11,...
%                                                   CxVals,CyVals,insideCutCell,ExtForcing,t0);
%     S_NL_ISO_debug(kj) = S_iso;
%     S_cutVals(kj) = S_cut(kj);
%     S_intVals(kj) = S_int(kj);
% 
% end
% diff_k = S_NL_ISO_debug - (S_intVals + S_cutVals);
% [maxAbsDiff, idx] = max(abs(diff_k));
% fprintf('max |S_iso - (S_int+S_cut)| = %e at kj=%d, k=%g\n', ...
%         maxAbsDiff, idx, kVals(idx));
% 
% G_debug = trapz(kVals, S_NL_ISO_debug);
% G_int = trapz(kVals, S_intVals);
% G_cut = trapz(kVals, S_cutVals);
% 
% fprintf('Interior-only sum ∑TΔk = %.6e\n', G_int);
% fprintf('Cut-cell-only sum ∑TΔk = %.6e\n', G_cut);
% fprintf('Total (check): %.6e\n', G_int + G_cut);
% fprintf('Inviscid, θ=1, simple E(k): ∑ T(k)Δk = %.6e\n', G_debug);
Eic = icGenerate(kVals);
loglog(kVals,Eic);
mu = getMu(Eic,kVals);
mu1 = A1*mu;
[weight,ySquaredAvg,triadFlag,outsideCutCell,insideCutCell,CxVals,CyVals,Q11]=midpoint2dShoelace(kVals);
kmin = 0.01;
[weight, CxVals,CyVals, dv, edges] = buildTriadWeightsCentroidsExact(kVals, kmin);
T = zeros([kLength,kLength,kLength]);
Ssym = T;
Tkj = T;
res = 0;
% for kj = 1:kLength
%     for pj = 1:kLength
%         for qj = 1:kLength
% triadFlag(pj,qj,kj) = (weight(pj,qj,kj) > 0);
% 
% 
%         end
%     end
% end
for kj = 1:length(kVals)
            kjl = kj;
            k = kVals(kjl);
          %[forcing,res,percentRes] = TwoD_MidpointTest(bf,nu,D,kjl,kVals(kjl),kVals,HPOL,HDIR,HT,E,F,ET,mu1,mu3,t,weight,du,qjMin,qjMax,ySquaredAvg,triadFlag,outsideCutCell,Q11,CxVals,CyVals,insideCutCell,ExtForcing,t0);
      % [S_NL_ISO,S_NL_DIR,S_NL_POL,ST_NL_ISO,ST_NL_DIR,SF_NL,phi33,Tkj,S_inter] = getTheIndividualTerms(0,nu,D,kj,k,kVals,HPOL,HDIR,HT,Eic,F,ET,mu1,mu3,tIC,weight,triadFlag,outsideCutCell,Q11,CxVals,CyVals,insideCutCell,ExtForcing,t0);
[S_NL_ISO,res2,Tkj]=conservativeS(kVals, edges, kj, E,CxVals,CyVals, weight);

    % forcing(1) = S_NL_ISO+0*N*F(kj)+ExtForcing(kj) ;
   res = max(res,abs(res2));
    S_NL_ISOicvals(kj) = S_NL_ISO;

   % S_interior(kj) = S_inter;
   T(kj,:,:) = Tkj(kj,:,:);
    
   
   
end
maxResidual = 0;
%symmetry check
% for kj = 2:kLength
%     for pj = 2:kLength
%         for qj = 2:kLength
% ok = triadFlag(pj,qj,kj)==1 && outsideCutCell(pj,qj,kj)==0 && insideCutCell(pj,qj,kj)==0 && ...
%      triadFlag(qj,kj,pj)==1 && outsideCutCell(qj,kj,pj)==0 && insideCutCell(qj,kj,pj)==0 && ...
%      triadFlag(kj,pj,qj)==1 && outsideCutCell(kj,pj,qj)==0 && insideCutCell(kj,pj,qj)==0;
% 
% if ok
% 
% 
%             %%residual = (T(kj,pj,qj) + T(pj,kj,qj) + T(qj,pj,kj))/(abs(min([dTk,dTp,dTq]))+1e-13);
%             %residual = T(kj,pj,qj) + T(pj,kj,qj);
%           residual = T(kj,pj,qj) + T(pj,qj,kj) + T(qj,kj,pj);
%             if abs(residual) > .0001
%                 fprintf('kj=%i,pj=%i,qj=%i \n',kj,pj,qj)
% 
%             end
%             maxResidualTemp = max([maxResidual,abs(residual)]);
%             if maxResidualTemp > maxResidual
% 
%             maxResidual = maxResidualTemp;
%             maxResidualLoc = [kj,pj,qj];
%             end
% end
%         end
%     end
% end
% maxResidual
% maxResidualLoc

% G_interior = trapz(kVals, S_interior);
% G_all      = trapz(kVals, S_NL_ISOicvals);



% E0 = Eic./(4*pi*kVals.^2);
% for kj = 1:kLength
%     for pj = 1:kLength
%         for qj = 1:kLength
%            T(kj,pj,qj) = E0(qj)*(E0(pj)-E0(kj));
%         end
%     end
% end


for kj = 1:kLength
    for pj = 1:kLength
        for qj = 1:kLength
            if triadFlag(qj,pj,kj) && triadFlag(pj,qj,kj)
               
            %Ssym(kj,pj,qj) = 0.5*(T(kj,pj,qj)+T(kj,qj,pj));
                Ssym(kj,pj,qj) = T(kj,pj,qj);
            end
        end
    end
end

resT=zeros([kLength kLength kLength]);
resTmax = 0.0;

for kj = 1:kLength
    for pj = 1:kLength
        for qj = 1:kLength
          if triadFlag(pj,qj,kj) && triadFlag(kj,pj,qj) && triadFlag(qj,kj,pj)

                        resT(kj,pj,qj) = Ssym(kj,pj,qj)+Ssym(pj,qj,kj)+Ssym(qj,kj,pj);
                        if abs(resT(kj,pj,qj)) > (abs(resTmax))
                            resTmax = resT(kj,pj,qj);
                            resTmaxCoords = [kj,pj,qj];
                        end

                    
                
           end

        end
    end
end
resTmax
resTmaxCoords
kj = resTmaxCoords(1); pj = resTmaxCoords(2); qj = resTmaxCoords(3);
resTmaxCoords
Cx1 = CxVals(pj,qj,kj) 
Cx2 = CxVals(qj,kj,pj) 
Cx3 = CxVals(kj,pj,qj)

Cy1 = CyVals(pj,qj,kj)
Cy2 = CyVals(qj,kj,pj)
Cy3 = CyVals(kj,pj,qj)
fprintf('Ssym(kj,pj,qj)=%f Ssym(pj,qj,kj)=%f Ssym(qj,kj,pj)=%f \n',Ssym(kj,pj,qj),Ssym(pj,qj,kj),Ssym(qj,kj,pj))
 figure()
 T_E_ICzero = trapz(kVals,S_NL_ISOicvals);
 T_E_ICpartialsTrapz = cumtrapz(kVals,(kVals.^2).*S_NL_ISOicvals);
 T_E_ICpartialsMid = midpoint1d(kVals,S_NL_ISOicvals)*0;
 dissipationIC = 2*nu*trapz(kVals,kVals.^2.0 .* Eic);
 dissipationICpartials = 2*nu*cumtrapz(kVals,kVals.^2.0 .* Eic);
 semilogx(kVals,T_E_ICpartialsTrapz,kVals,T_E_ICpartialsMid,kVals,-dissipationICpartials)
 legend(['T_{E,iso} partial sums at IC trapz     '; ...
         'T_{E,iso} partial sums at IC midpoint  ';
         'KEdissipation partial sums             '])
 title('Initial condition partial sums')
 
figure()
 semilogx(kVals, S_NL_ISOicvals+ExtForcing-2*nu*kVals.^2 .*Eic)
 title('IC E ISO balance')
 axis([-inf inf -inf inf ])

 symErr = max(abs(weight - permute(weight,[2 1 3])),[],'all');
fprintf('symErr=%g\n', symErr);

A_k = zeros(length(kVals),1);
for kj=1:length(kVals)
    A_k(kj) = sum(weight(:,:,kj),'all');
end
fprintf('A_k min/max = %g  %g\n', min(A_k), max(A_k));

TE = 4*pi*kVals(:).^2 .* S_NL_ISOicvals(:);
fprintf('trapz(TE)=%g, trapz(S_NL_ISO)=%g\n', trapz(kVals,TE), trapz(kVals,S_NL_ISOicvals));

[S_NL_E, diag] = transfer_scatter_add_kernelE0(kVals, edges, E, weight, CxVals, CyVals);
fprintf('FV total energy transfer = %.20e\n', diag.FV_total_energy_transfer);
fprintf('max triad energy residual = %.20e\n', diag.maxTriadEnergyResidual);
fprintf('triads used = %d\n', diag.numTriadsUsed);

partial = cumsum(S_NL_E .* diff(edges));
semilogx(sqrt(edges(1:end-1).*edges(2:end)), partial), grid on, yline(0,'--')
title('FV partial sums of energy transfer (scatter-add, your kernel)')


dk   = diff(edges(:));
krep = sqrt(edges(1:end-1).*edges(2:end));   % log-midpoint per bin

TEbin = 4*pi*krep.^2 .* S_NL_ISOicvals(:);   % if S_NL_ISO is in E0-units
Total = sum(TEbin .* dk);
figure()
fprintf('FV total energy transfer = %.6e\n', Total);
partial = cumsum(TEbin .* dk);
semilogx(krep, partial), grid on, yline(0,'--')
title('FV partial sums of energy transfer')
a = -1;

bad_qk = 0; bad_kp = 0; tot = 0;
Cx = CxVals;
Cy = CyVals;

for i = 2:kLength-1
  for j = 2:kLength-1
    for l = 2:kLength-1
      % only test where all cyclic perms exist
      if ~( triadFlag(i,j,l) && triadFlag(j,l,i) && triadFlag(l,i,j) )
        continue
      end
      tot = tot+1;

      % pq|k
      pstar = Cx(i,j,l); qstar = Cy(i,j,l); kfix = kVals(l);
      ok_pq = (abs(pstar-qstar) <= kfix+1e-12) && (kfix <= pstar+qstar+1e-12);

      % qk|p
      qstar2 = Cx(j,l,i); kstar2 = Cy(j,l,i); pfix = kVals(i);
      ok_qk = (abs(qstar2-kstar2) <= pfix+1e-12) && (pfix <= qstar2+kstar2+1e-12);

      % kp|q
      kstar3 = Cx(l,i,j); pstar3 = Cy(l,i,j); qfix = kVals(j);
      ok_kp = (abs(kstar3-pstar3) <= qfix+1e-12) && (qfix <= kstar3+pstar3+1e-12);

      if ~ok_pq
        error('pq check failed: your Cx,Cy do not even satisfy pq|k');
      end
      if ~ok_qk, bad_qk = bad_qk+1; end
      if ~ok_kp, bad_kp = bad_kp+1; end
    end
  end
end

fprintf('tested=%d, bad_qk=%d (%.2f%%), bad_kp=%d (%.2f%%)\n', ...
        tot, bad_qk, 100*bad_qk/tot, bad_kp, 100*bad_kp/tot);



mu= getMu(E,kVals); %mu3 of 'Spectral modelling for passive scalar dynamics in hom anis turb'
            mu1 = A1*mu;
            mu3 = A3*mu;
[HDIR,HPOL,HT] = getHvals(kLength,E,ET,EHdir,EHpol,ETH);
t0 = 50.0*TauL0
for kj = 1:length(kVals)
            kjl = kj;
            k = kVals(kjl);
          %[forcing,res,percentRes] = TwoD_MidpointTest(bf,nu,D,kjl,kVals(kjl),kVals,HPOL,HDIR,HT,E,F,ET,mu1,mu3,t,weight,du,qjMin,qjMax,ySquaredAvg,triadFlag,outsideCutCell,Q11,CxVals,CyVals,insideCutCell,ExtForcing,t0);
       [S_NL_ISO,S_NL_DIR,S_NL_POL,ST_NL_ISO,ST_NL_DIR,SF_NL,phi33] = getTheIndividualTerms(bf,nu,D,kj,k,kVals,HPOL,HDIR,HT,E,F,ET,mu1,mu3,t,weight,triadFlag,outsideCutCell,Q11,CxVals,CyVals,insideCutCell,ExtForcing,t0);

    % forcing(1) = S_NL_ISO+a*N*F(kj)+ExtForcing(kj) ;
   
    ST_NL_ISOvals(kj) = ST_NL_ISO;
    Ftimes2bf(kj)= 2.*bf*F(kj);
   
   
end
a = 1;
kVals = [0 kVals(a:end)];
ETnew = [0 ETnew(a:end)];
ST_NL_ISOvals=[0 ST_NL_ISOvals((a:end))];
Ftimes2bf = [0 Ftimes2bf(a:end)];
F = [0 F(a:end)];
kaiPartials = 2.0*D*cumtrapz(kVals,kVals.^2.*ETnew);
ST_NL_ISOpartials = cumtrapz(kVals,ST_NL_ISOvals);
Ftimes2bfPartials = cumtrapz(kVals,Ftimes2bf);
    
figure();
semilogx(kVals,ST_NL_ISOpartials,'--g')
hold on;
semilogx(kVals,Ftimes2bfPartials,'w')
hold on;
semilogx(kVals,-kaiPartials)
legend(['T_{NL} isotropic NL transfer '; ...
        '2*N*buoyancyFlux             '; ...
        '-2*\chi                      '],'Location','east')
title('logarithmic x axis plot of partial integral sums (\int 0 to k)')
xlabel('log k')
%===================================================================================
figure();

loglog(kVals,-ST_NL_ISOpartials,'--g')
hold on;
loglog(kVals,Ftimes2bfPartials,'w')
hold on;
loglog(kVals,kaiPartials)
legend(['-T_{NL} isotropic NL transfer'; ...
        'buoyancyFlux 2*bf*F          '; ...
        '2*\chi                       '],'Location','southeast')
title('loglog plot of partial integral sums (\int 0 to k)')
xlabel('log k')

figure()
loglog(kVals,abs(ST_NL_ISOpartials-kaiPartials+Ftimes2bfPartials+[0 ExtForcing]))
title('T_{NL} - 2\chi+ 2*bf*F partial integral sums')
xlabel('k')

figure()
loglog(kVals,2*D*kVals.^2.*ETnew,kVals,2*bf*F,'w',kVals,-ST_NL_ISOvals,'--g')
title('loglog plot of 2Dk^2ET(k),2bf*F(k), -ST_{NL}^{ISO}(k)')
legend(['2*D*k^2*ET(k)'; ...
        'Bflux        '; ...
        'T_{NL}       '],'Location','southeast')
axis([-inf inf 10^-10 inf])
figure()
semilogx(kVals,2*D*kVals.^2.*ETnew,kVals,2*bf*F,'w',kVals,ST_NL_ISOvals,'--g')
title('log x axis plot of 2Dk^2ET(k),2bf*F(k), -ST_{NL}^{ISO}(k)')
legend(['2*D*k^2*ET(k)'; ...
        'Bflux        '; ...
        'T_{NL}       '],'Location','southeast')
figure()

semilogx(kVals,Ftimes2bf)

title('logx plot of 2*bf*F')


L1 = 0; L2 = 0;
L3 = 0;
for i = 1:kLength
HPOLmatrix = [-HPOL(i)/2 0 0; 0 -HPOL(i)/2 0; 0 0 HPOL(i)];
L1 = max(L1,abs(eigs(HPOLmatrix,1)));
L1vals(i)=abs(eigs(HPOLmatrix,1));
HDIRmatrix = [-HDIR(i)/2 0 0; 0 -HDIR(i)/2 0; 0 0 HDIR(i)];
L2 =  max(L2,abs(eigs(HDIRmatrix,1)));
L2vals(i)=abs(eigs(HDIRmatrix,1));
HTmatrix = [-HT(i)/2 0 0; 0 -HT(i)/2 0; 0 0 HT(i)];
L3 =  max(L3,abs(eigs(HTmatrix,1)));
L3vals(i) = abs(eigs(HTmatrix,1));



end
L1
L2
L3
figure()
semilogx(kVals(2:end),L1vals,kVals(2:end),L2vals,kVals(2:end),L3vals);
legend(['sup_{i} (Eig_{i} (H^{Poloidal})    ';
        'sup_{i} (Eig_{i} (H^{Directional}) ';
        'sup_{i} (Eig_{i} (H_{T}^{Poloidal})'])
title('anisotropy eigenvalues of spectras vs logk')
axis([-inf inf -inf 0.5])
end




