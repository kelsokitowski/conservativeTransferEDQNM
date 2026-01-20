function dA = cutBykPlusP(p,q,k,pM,pP,qM,qP)

%this function returns area that includes center point.
%q = k+p or q-k-p = 0 or p = q-k
c1s = (qP-k-pM); %corner 1 sign (upper left, numbered clockwise)
                    c2s = (qP-k-pP);
                    c3s = (qM-k-pP);
                    c4s = (qM-k-pM);
                    centerSign = q-k-p;
                    
                    s1Flag = 0;
                    s2Flag = 0;
                    s3Flag = 0;
                    s4Flag = 0;
                    x = zeros([1,8]);
                    y = x;
                    if c1s*c2s<0 %then intersection between corners
                        star1 = [qP-k,qP];
                        s1Flag = 1;
                    end
                    if c2s*c3s<0
                        star2 = [pP,k+pP];
                        s2Flag = 1;
                    end
                    if c3s*c4s <0
                        star3 = [qM-k,qM];
                        s3Flag = 1;
                    end
                    if c4s*c1s<0
                        star4 = [pM,k+pM];
                        s4Flag = 1;
                    end
                    index = 0;
                    if c1s*centerSign > 0 
                        index = index+1;
                        y(index) = qP;
                        x(index) = pM;
                    end
                    if (s1Flag == 1)
                    index = index + 1;
                    x(index) = star1(1);
                    y(index) = star1(2);
                    end
                    if c2s*centerSign > 0 
                        index = index+1;
                        y(index) = qP;
                        x(index) = pP;
                    end
                    if (s2Flag == 1)
                    index = index + 1;
                    x(index) = star2(1);
                    y(index) = star2(2);
                    end
                    if c3s*centerSign > 0 
                        index = index+1;
                        y(index) = qM;
                        x(index) = pP;
                    end
                    if (s3Flag == 1)
                    index = index + 1;
                    x(index) = star3(1);
                    y(index) = star3(2);
                    end
                    if c4s*centerSign > 0 
                        index = index+1;
                        y(index) = qM;
                        x(index) = pM;
                    end
                     if (s4Flag == 1)
                    index = index + 1;
                    x(index) = star4(1);
                    y(index) = star4(2);
                     end
                dA = polyarea(x(1:index),y(1:index));
end