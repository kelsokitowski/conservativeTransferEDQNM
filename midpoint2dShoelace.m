function [areas,ySquaredAvg,triadFlag,outsideCutCell,insideCutCell,CxVals,CyVals,Q11]=midpoint2dShoelace(kVals)
%midpoints = kVals;



%k= kVals(80);
kLength = length(kVals);
areas = zeros([kLength kLength kLength]);
areas2=areas;
triadFlag = areas;
ySquaredAvg = areas;
xAvg = areas;
yAvg = areas;
zAvg = areas;
totallyInside = areas;
outsideCutCell = areas;
insideCutCell = areas;
Q11 = zeros([2 kLength kLength kLength]);
CyVals = areas;
CxVals = areas;
for kj = 1:kLength
    k=kVals(kj);
    triadWeights = zeros([kLength kLength]); %for visualization only
    tic
    for qj = 1:kLength
        
        for pj = qj:kLength
            q = kVals(qj);
            p = kVals(pj);
          
               [pM,pP,qM,qP] =getCorners(pj,qj,kLength,p,q,kVals);
                A_rectangle = (pP-pM)*(qP-qM);
                pCoords = [pM,pP];
                qCoords = [qM,qP];
                dAcutByKPlusP = 0.0;
                dAcutByKMinusP = 0.0;
                dAcutByPMinusK = 0.0;
                anyCornerInside = 0;
                for i = 1:2
                    for j = 1:2
                        q0 = qCoords(i);
                        p0 = pCoords(j);
                        anyCornerInside = anyCornerInside + triadCondition(q0,p0,k);

                    end
                end
                if anyCornerInside > 0 | triadCondition(q,p,k)
                    if kj == 1 && pj ==92 && qj == 91
%    disp('debug here')
   %% hold on; scatter3(pj,qj,15,'yellow')
end
                    [Cx,Cy,A] =centroidOfEdgeCell(qM,qP,pM,pP,k);
                    CxVals(pj,qj,kj) = Cx;
                    CyVals(pj,qj,kj) = Cy;
                    Q11 = boundingBoxInterpolation(Cx,Cy,kj,pj,qj,p,q,Q11);
                    areas2(pj,qj,kj) = A;
		    if A<0.00000000001
			    disp('areas2 somehow zero :58 of midpoint2dShoelace.m')
		    end
                    if anyCornerInside == 4
                        %[Cx,Cy,A] =centroidOfEdgeCell(qM,qP,pM,pP,k);
                       testCheck =  A-A_rectangle;
                       if abs(testCheck) > 1e-12
                           disp('inconsistency in centroidOfEdgeCell area')
                       end
                    end
                end
                flag1 = 0;
                flag2 = 0;
                flag3 = 0;
                centerSign = q-p-k; %q = k+p
                for i = 1:2
                    for j = 1:2


                        p0 = pCoords(i);
                        q0 = qCoords(j);
                        sign = q0-p0-k;
                        if (centerSign*sign <= 0) %upper domain intersection
                           
                            dAcutByKPlusP = cutBykPlusP(p,q,k,pM,pP,qM,qP);
                            if dAcutByKPlusP > A_rectangle
                                fprintf('Cut area exceeds rectangle at pj=%d,qj=%d,kj=%d\n', pj,qj,kj);
                            end
                            dAcutByKPlusP = A_rectangle-dAcutByKPlusP;
                            if dAcutByKPlusP < 0
                                fprintf('NEGATIVE intermediate value at pj=%d,qj=%d,kj=%d\n', pj,qj,kj);
                            end
                          
                            flag1 = 1;
                        end


                    end
                end



                centerSign = q-k+p;  %q=k-p
                for i = 1:2
                    for j = 1:2


                        p0 = pCoords(i);
                        q0 = qCoords(j);
                        sign = q0-k+p0;
                        if (centerSign*sign <= 0) %lower left domain intersection
                            dAcutByKMinusP = cutBykMinusP(p,q,k,pM,pP,qM,qP);
                            dAcutByKMinusP = A_rectangle-dAcutByKMinusP;
                            flag2 = 1;
                        end


                    end
                end




                centerSign = q+k-p;  %q=p-k
                for i = 1:2
                    for j = 1:2


                        p0 = pCoords(i);
                        q0 = qCoords(j);
                        sign = q0+k-p0;
                        if (centerSign*sign <= 0) %lower right domain intersection
                            dAcutByPMinusK = cutByPMinusK(p,q,k,pM,pP,qM,qP);
                            dAcutByPMinusK = A_rectangle-dAcutByPMinusK;
                            flag3 = 1;
                        end


                    end
                end
                if triadCondition(q,p,k)
                    areas(pj,qj,kj)= A_rectangle-flag1* dAcutByKPlusP-flag2*dAcutByKMinusP-flag3*dAcutByPMinusK;
                    triadFlag(pj,qj,kj) = 1;
                    
                    if (flag1+flag2+flag3 < 0.5) %inside and no cut
                        totallyInside(pj,qj,kj) = 1;
                    else
                        insideCutCell(pj,qj,kj) = 1;
                    end
                elseif (flag1+flag2+flag3 > 0.5) %then outside but piece inside, so complement of above is the areas inside
                    areas(pj,qj,kj) = flag1* dAcutByKPlusP+flag2*dAcutByKMinusP+flag3*dAcutByPMinusK;
                    triadFlag(pj,qj,kj) = 1;
                    outsideCutCell(pj,qj,kj) = 1;
                    
                end
              
                % xAvg(qj,pj,kj) = xAvg(pj,qj,kj) ;
                % yAvg(qj,pj,kj) = yAvg(pj,qj,kj) ;
                % zAvg(qj,pj,kj) = zAvg(pj,qj,kj) ;
                % ySquaredAvg(qj,pj,kj) = ySquaredAvg(pj,qj,kj) ;
               if abs(abs(areas2(pj,qj,kj))-abs(areas(pj,qj,kj))) >1.0e-04
                   disp('uhoh areas2(pj,qj,kj) neq areas(pj,qj,kj) :144 midpoint2dShoelace')
		  fprintf('pj=%i,qj=%i,kj=%i ',pj,qj,kj)
          fprintf('flag1=%i flag2=%i flag3=%i \n',flag1,flag2,flag3)
		  fprintf('areas2(pj,qj,kj)=%f,areas(pj,qj,kj)=%f',areas2(pj,qj,kj),areas(pj,qj,kj))
          a2MinusA = abs(abs(areas2(pj,qj,kj))-abs(areas(pj,qj,kj)));	
		   fprintf('outsideCutCell(pj,qj,kj)=%i, a2-A=%f \n',outsideCutCell(pj,qj,kj),a2MinusA)
			
               end
                areas(qj,pj,kj) = areas(pj,qj,kj) ; %by 45 degree symmetry
                triadFlag(qj,pj,kj) = triadFlag(pj,qj,kj) ;
                totallyInside(qj,pj,kj) = totallyInside(pj,qj,kj);
                outsideCutCell(qj,pj,kj) = outsideCutCell(pj,qj,kj);
                insideCutCell(qj,pj,kj) =  insideCutCell(pj,qj,kj);
                areas2(qj,pj,kj) = areas2(pj,qj,kj);
                


                
                

              

            if triadCondition(q,p,k)
                triadWeights(pj,qj)= 1; %for domain visualization only

            end

        end

    end
    if kj == 1
%visualizeDomain(kVals,triadWeights,0,0,kj)
    
    end
    %end



    %GET CENTROIDS AND INTERPOLATION COORDS

    for qj = 1:kLength

        for pj = 1:kLength
            k = kVals(kj);
            q = kVals(qj);
            p = kVals(pj);
 if kj == 133 && qj == 133 && pj == 1
 %    disp('debug')
 end
            [pM,pP,qM,qP] =getCorners(pj,qj,kLength,p,q,kVals);
            A_rectangle = (pP-pM)*(qP-qM);
            pCoords = [pM,pP];
            qCoords = [qM,qP];
            anyCornerInside = 0;
            for i = 1:2
                for j = 1:2
                    q0 = qCoords(i);
                    p0 = pCoords(j);
                    anyCornerInside = anyCornerInside + triadCondition(q0,p0,k);

                end
            end
            
            if anyCornerInside > 0 || triadCondition(q,p,k)
                %kj
                %pj
                %qj
                if kj == 1 && pj == 132 && qj == 133
%    disp('debug here')
end
                [Cx,Cy,A] =centroidOfEdgeCell(qM,qP,pM,pP,k);
                CxVals(pj,qj,kj) = Cx;
                CyVals(pj,qj,kj) = Cy;
                Q11 = boundingBoxInterpolation(Cx,Cy,kj,pj,qj,p,q,Q11);
                % areas2(pj,qj,kj) = A;
                if anyCornerInside == 4
                    [Cx,Cy,A] =centroidOfEdgeCell(qM,qP,pM,pP,k);
  %                  testCheck =  A-A_rectangle
                end
            end
        end
    end


    kj
    toc
end

% for kj = 1:kLength
%     for qj = 1:kLength
%         for pj = 1:kLength
%             % CxVals(qj,pj,kj) = CxVals(pj,qj,kj);
%             % CyVals(qj,pj,kj) = CyVals(pj,qj,kj);
%              CxValsAvg(pj,qj,kj) = ( CxVals(qj,pj,kj) + CxVals(pj,qj,kj) )/2;
%              CyValsAvg(pj,qj,kj) = ( CyVals(qj,pj,kj) + CyVals(pj,qj,kj) )/2;
% 
% 
%         end
%     end
% end
% 
% CxVals = CxValsAvg;
% CyVals = CyValsAvg;


% Cx = CxVals(:,:,kj);
% Cy = CyVals(:,:,kj);
% 
% CxT = Cy.';   % what Cx should be under swap
% CyT = Cx.';   % what Cy should be under swap
% 
% CxVals(:,:,kj) = 0.5*(Cx + CxT);
% CyVals(:,:,kj) = 0.5*(Cy + CyT);

%if 1 == 0
for kj = 1:kLength

    %debug test vs analytical solution
    integral = 0;
    integralAreas2 = 0;
    integralExcludingOutsideCenters = 0.0;
    integralExcludingOutsideCenters2 = 0.0;
    %kj
    k = kVals(kj);
    for pj=1:kLength
        for qj = 1:kLength
            p = kVals(pj);
            q = kVals(qj);
           
                integral = integral + areas(pj,qj,kj)*1;
                integralAreas2 = integralAreas2 + areas2(pj,qj,kj)*1;
                if triadCondition(q,p,k) %check performance when excluding cell centers outside domain
                integralExcludingOutsideCenters =  integralExcludingOutsideCenters + ...
                    areas(pj,qj,kj)*1;
                integralExcludingOutsideCenters2 =  integralExcludingOutsideCenters2 + ...
                    areas2(pj,qj,kj)*1;
                end




                if (pj==6 && qj == 1 && kj == 1)
   %                 disp('check here')
                end
                %\\\\\\\\COMPUTE AVERAGE TRIG QUANTITIES
                % if 1==0
                % %if triadFlag(pj,qj,kj) == 1
                % if outsideCutCell(pj,qj,kj)==1
                %     %now compute average cosines x,y,z over domain
                %     if (totallyInside(pj,qj,kj) < 0.5) %then at least one corner of cell outside domain
                %         [pM_cell,pP_cell,qM_cell,qP_cell] =getCorners(pj,qj,kLength,p,q,kVals);
                % 
                %         % nPPpoints = min(ceil(max((pP_cell-pM_cell)/(kVals(2)-kVals(1)),4)),calcBound)+1;
                %         nPPpoints = max(ceil((pP_cell-pM_cell)/(kVals(2)-kVals(1)))*4,1000)+1;
                % 
                % 
                %         pRange = linspace(pM_cell,pP_cell, nPPpoints);
                % 
                % 
                %         for pb = 1:length(pRange)
                %             pLocal = pRange(pb);
                %             outFlag = 1; %if 1, portion of vertical slice in domain
                % 
                %             if triadCondition(qM_cell,pRange(pb),k)
                % 
                %                 qM = qM_cell;
                %                 qP = min(qP_cell,pLocal+k); %in case top boundary cuts here.
                %             else
                %                 if ( qM_cell < k-pLocal ) %then we are below left bottom boundary
                %                     qM = k-pLocal;
                %                     if qM>=qP_cell  %then there is no vertical slice in domain here.
                %                         outFlag = 0; %do not do the dp integral.
                %                     else
                %                         qP = min(qP_cell,pLocal+k); %in case top boundary cuts here.
                %                     end
                % 
                %                 elseif (qM_cell < pLocal-k)
                %                     qM = pLocal - k;
                %                     if qM >= qP_cell %then there is no vertical slice in domain here.
                %                         outFlag = 0;
                %                     else
                %                         qP = min(qP_cell,pLocal+k); %in case top boundary cuts here.
                %                     end
                %                 elseif (qM_cell>k+pLocal) %bottom of vertical slice is already above domain.
                %                     outFlag = 0; %0*dp here for final integral.
                %                 else
                %                     disp('fallthrough, uhoh')
                %                 end
                %             end
                % 
                % 
                % 
                %             if outFlag == 0
                %                 kernel1(pb) = 0.0;
                %             else
                %                 kernel1(pb) = -(3*k^4*qM - 3*k^4*qP - 6*k^2*pLocal^2*qM + 6*k^2*pLocal^2*qP + 6*k^2*qM^2*qP - 6*k^2*qM*qP^2 + 3*pLocal^4*qM - 3*pLocal^4*qP - 6*pLocal^2*qM^2*qP + 6*pLocal^2*qM*qP^2 + qM^4*qP - qM*qP^4)/(12*k^2*qM*qP);
                %             end
                %             %above is the analytically integrated y^2 along q
                %             %(vertical),which we need to determine qM,qP for a
                %             %given p.
                % 
                % 
                % 
                %         end
                %         ySquaredKernel = trapezoidalIntegration(pRange,kernel1(1:length(pRange)));
                % 
                % 
                %     elseif (totallyInside(pj,qj,kj) > 0.5)  %totallyInside, use analytical expression
                % 
                %         [pM,pP,qM,qP] =getCorners(pj,qj,kLength,p,q,kVals);
                %         ySquaredKernel = -((3*qM - 3*qP)*(-pM^5 + pP^5))/(60*k^2*qM*qP) - ((-6*k^2*qM + 6*k^2*qP - 6*qM^2*qP + 6*qM*qP^2)*(-pM^3 + pP^3))/(36*k^2*qM*qP) - ((3*k^4*qM - 3*k^4*qP + 6*k^2*qM^2*qP - 6*k^2*qM*qP^2 + qM^4*qP - qM*qP^4)*(pP - pM))/(12*k^2*qM*qP);
                % 
                %     end
                % 
                %     ySquaredAvg(pj,qj,kj) = ySquaredKernel/areas(pj,qj,kj);
                %     if 1-ySquaredAvg(pj,qj,kj)<0
                %         if abs(1-ySquaredAvg(pj,qj,kj)) < 1e-9
                %             ySquaredAvg(pj,qj,kj)=1;
                %         else
                %             disp('negative value at')
                %             pj
                %             qj
                % 
                % 
                %             visualizeDomain(kVals,triadWeights,1,[kVals(pj),kVals(qj)],kj)
                %             ySquaredAvg(pj,qj,kj) = 1;
                %         end
                %     end
                % end  
                % end %if (1==0)  (use to turn on or off averaging.)
                %\\\\\ END OF COMPUTING AVERAGE TRIG QUANTITIES

          
        end
    end
    analytical = 2*k*max(kVals)-1.5*k^2;

    res(kj) = integral-analytical;
    res2(kj) = integralAreas2-analytical;
    percentError(kj) = abs(res(kj)/analytical*100);
    resNoOutside(kj)=integralExcludingOutsideCenters-analytical;
    resNoOutside2(kj)=integralExcludingOutsideCenters2-analytical;
    percentErrorNoOutside(kj) = abs(resNoOutside(kj)/analytical*100);
    percentErrorNoOutside2(kj) = abs(resNoOutside2(kj)/analytical*100);
    
    %visualizeDomain(kVals,triadWeights,0,0,kj)
  


end
% load("weight141bisector.mat")
% find(abs(weight(:,:,:)-areas(:,:,:))>.00000001)




%areas = weight;

maxPercentError = max(percentError)
maxRes = max(res)

maxPercentErrorNoOutside = max(percentErrorNoOutside)
maxResNoOutside = max(resNoOutside)
figure()
subplot(1,3,1)
loglog(kVals,abs(res))
title('Residual of Integral of 1*dA')
xlabel('k')

subplot(1,3,2)
loglog(kVals,percentErrorNoOutside,'xw',kVals,percentErrorNoOutside2,'--g')
title('Percent Error of Integral of 1*dA with removal of outside centers')
xlabel('k')

subplot(1,3,3)
loglog(kVals,abs(percentError))
title('Percent error of Integral of 1*dA')
xlabel('k')


% disp('end of function midPoint2dShoelace')

%check bilinear stencil p,q symmetry
for kj = 1:kLength
    for pj = 1:kLength
        for qj = 1:kLength
            if (Q11(1,pj,qj,kj)~= Q11(2,qj,pj,kj))
                disp('symmetry error')
            end
            if (Q11(2,pj,qj,kj)~= Q11(1,qj,pj,kj))
                disp('symmetry error')
            end
        end
    end
end








end





%remove the areas inside associated with points outside of the domain
%(although it shouldn't be called there anyways..)
% for i  = 1:kLength
%     for j = 1:kLength
%         for m = 1:kLength
%             k = kVals(i);
%             p = kVals(j);
%             q = kVals(m);
%             if (~triadCondition(q,p,k))
%                 areas(i,j,m)=0.0;
%             end
% 
%         end
%     end
% end
% end
