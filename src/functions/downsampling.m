function [M_kernelR, indr, intind] = downsampling(nx, maxn, int_tol, M_kernel, ts_index, k)
    M_kernelR = cell(1, nx);
    intind = cell(1,nx);
    indr = cell(1,nx);
    lr_index = zeros(1,length(k));
    for i_xx = 1:nx
              
               ref = cumsum(M_kernel{i_xx}(:,1)); %118-119 redundant
               
    %           testdiff = diff([0;inddex]);       % add a zero on the top of inddex and get their neighbouring elements difference
    %           inddex = inddex - (testdiff - 1);
               ishift = 1;
               subind = ishift;
               maxnT = maxn;
               % decide if the approximating range being too large, if maxnT =
               % 1, no reducing action will be made
                if ishift + (2*maxn+1) > ts_index(i_xx)
                    disp('maxn is too big for this k, setting to 1')
                    maxnT = 1;
                end
               
               startnT = 1;
               
               while ishift + (2*maxnT+1) <= ts_index(i_xx)
                    if startnT < maxnT
                        for n = startnT:maxnT
                            subindT = [subind;ishift+2*n+1]; % n index -> odd -> 
                            %inddexT = inddex(subindT);
                            testdiffT = diff([0;subindT]);
                            inddexT = subindT - (testdiffT - 1)/2;
                            
                            interr = abs(sum(testdiffT.*M_kernel{i_xx}(inddexT,1)) ...
                                 - ref(2*n+1 + ishift))./ref(2*n+1 + ishift); 
                            if n == maxnT && interr < 0.25*int_tol
                                interr = 0;  
                            end
                            
                             
                            if interr > int_tol % odd number step
                                break           % the tolerance reached, such that jump out the approximation cycle
                            end
                        end
                        n = n - 1;
                        subind = [subind;ishift + 2*n+1];
                        %inddexT = inddex(subind);
                        testdiffT = diff([0;subind]);
                        inddexT = subind - (testdiffT - 1)/2; % indexing those ... back a few steps? corrector?
                        ishift = ishift + (2*n+1) ;
                        if ishift + (2*maxnT+1) > ts_index(i_xx)
                           maxnT = (- ishift + ts_index(i_xx) - 1)/2;
                           if maxnT < 1
                               break
                           end
                        end
                        startnT = n; %!!!!
                        if interr == 0 && startnT + 1 == maxnT
                           startnT = maxnT; 
                        end
                        if startnT > maxnT
                            startnT = floor(maxnT);
                        end
                    else
                        n = maxnT;
                        nf = floor(((ts_index(i_xx) - ishift)   -  1)/(2*n));  
                        nv = (n*ones(nf,1));
                        subind = [subind; ishift + cumsum(2*nv+1) ];
                        subind(subind > ts_index(i_xx)) = [];
                        %inddexT = inddex(subind);
                        testdiffT = diff([0;subind]); % an odd number, max = maxn*2 + 1
                        inddexT = subind - (testdiffT - 1)/2; % testdiffT - diff([0; inddexT]) = [0,..,0] almost
                        %inddexT_test = subind;
                        break
                    end
               end
         
    
               indr{i_xx}(:,1) = flip(inddexT);  % store reduced index with the cell-array form, from the last element to the first element - flipped
               intind{i_xx}(:,1) = testdiffT;      
               M_kernelR{i_xx}(:,1) = flip(intind{i_xx}(:,1).*M_kernel{i_xx}(inddexT,1)); %?, flipped and accumulated
               lr_index(1,i_xx) = length(inddexT);
    end
    disp('Compression') % compare kernel nodes; get 90% zero point off
    1 - mean(sum(lr_index)/(sum(ts_index)));
end
