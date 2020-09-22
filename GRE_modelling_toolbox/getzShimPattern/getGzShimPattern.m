function [ MzShim, GzShim] = getGzShimPattern( Gz_vec, te, pattern, start_echo)
%GETGZSHIMPATTERN Summary Returns the compensation moment of the z-Shim
%between ech echo based on the selected pattern. 
%   @Gz_vec         %Vector containing field gradients values (mT/m)
%   @pattern        %selected pattern


    Nte = length(te); 

    Gc = Gz_vec;
    Nslices = size(Gc,1); 
    switch pattern
       case '+A_-A_0_+B_-B'
            for i=1:Nslices
                %Check if Gz_low or Gz_high are zero 
                if sum(Gc(i,:)==0)  == 1
                    idx = find(Gc(i,:) ~= 0); 
                    for j=start_echo+1:2:Nte
                        
                        MzShim(i,j) =  Gc(i,idx).*te(j)*1E3;
                        MzShim(i,j+1) = - MzShim(i,j);
                    end

                else
                    
                    MzShim(i,:) = zeros(1,Nte); 
                    
                    MzShim(i,start_echo+1:4:end) = te(start_echo+1:4:end)*1E3*Gc(i,1);
                    MzShim(i,start_echo+2:4:end) = - MzShim(i,start_echo+1:4:end) ;
                    
                    MzShim(i,start_echo+3:4:end) = te(start_echo+3:4:end)*1E3*Gc(i,2);
                    MzShim(i,start_echo+4:4:end) = - MzShim(i,start_echo+3:4:end);

                end
            
            end
            
        case '+A_-A+B_-B_+A'
            
            for i=1:Nslices
                if sum(Gc(i,:)==0)  == 1
                    idx = find(Gc(i,:) ~= 0); 
                    for j=start_echo+2:Nte
                        if j== start_echo+2
                             MzShim(i,j) = Gc(i,idx).*te(start_echo+2)*1E3;
                        else
                            MzShim(i,j) = Gc(i,idx)*(te(j) -te(j-1))*1E3;
                        end

                    end

                else

                    for j=start_echo+1:2:Nte-1
                        dMa = te(j)*1E3*Gc(i,1); 

                        if j==start_echo+1
                            dMb = te(j+1)*1E3*Gc(i,2); 
                            MzShim(i,j) = dMa; 
                            MzShim(i,j+1) = -dMa + dMb; 
                        else
                            dMb = te(j-1)*1E3*Gc(i,2);      
                            MzShim(i,j) = dMa - dMb; 

                            dMb = te(j+1)*1E3*Gc(i,2);      
                            MzShim(i,j+1) = -dMa + dMb;  
                        end
                    end
                end
            end 
        
        case '+A/3_2A/3_A_-A_+B/3_+2B/3_B'
            
          for i=1:Nslices
               count = 1; 
               for j=start_echo+1:Nte
                   if count==1
                      MzShim(i,j) = 0.5*Gc(i,1).*te(j)*1E3;
                   elseif count==2 
                       MzShim(i,j) = 1*Gc(i,1).*te(j)*1E3 - MzShim(i,j-1);
                   elseif count==3
                       MzShim(i,j) = 1.5*Gc(i,1).*te(j)*1E3 -  MzShim(i,j-1) -  MzShim(i,j-2);
                   elseif count==4
                        MzShim(i,j) = -MzShim(i,j-1) -MzShim(i,j-2) -MzShim(i,j-3);
                   elseif count==5
                        MzShim(i,j) = 0.5*Gc(i,2).*te(j)*1E3;
                   elseif count==6 
                        MzShim(i,j) = 1*Gc(i,2).*te(j)*1E3 - MzShim(i,j-1);
                   elseif count==7
                        MzShim(i,j) = 1.5*Gc(i,2).*te(j)*1E3  -  MzShim(i,j-1) -  MzShim(i,j-2);
                   elseif count==8
                        MzShim(i,j) = -MzShim(i,j-1) -MzShim(i,j-2)- MzShim(i,j-3); 
                   end



                   if count ~= 8;
                     count = count + 1;
                   else
                     count = 1; 
                   end

               end
          end
          
          
        case  'Gc_split'

             for i=1:Nslices
               count = 1; 
               for j=start_echo+1:Nte
                   %check if a compensation gradient is 0
                   if sum(Gc(i,:)==0)  == 1
                       idx = find(Gc(i,:) ~= 0); 
                       if count==1
                          MzShim(i,j) = 0.5*Gc(i,idx).*te(j)*1E3;
                       elseif count==2 
                           MzShim(i,j) = 0.75*Gc(i,idx).*te(j)*1E3 - MzShim(i,j-1);
                       elseif count==3
                           MzShim(i,j) = 1*Gc(i,idx).*te(j)*1E3 -  MzShim(i,j-1) -  MzShim(i,j-2);
                       elseif count==4
                             MzShim(i,j) = 1.25*Gc(i,idx).*te(j)*1E3 -  MzShim(i,j-1) -  MzShim(i,j-2)- MzShim(i,j-3);
                       elseif count==5
                             MzShim(i,j) = 1.5*Gc(i,idx).*te(j)*1E3 -  MzShim(i,j-1) -  MzShim(i,j-2)- MzShim(i,j-3) - MzShim(i,j-4);
                       elseif count==6 
                             MzShim(i,j) = -MzShim(i,j-1) -MzShim(i,j-2)- MzShim(i,j-3) -  MzShim(i,j-4) - MzShim(i,j-5);  
                       end
                       
                       if count ~= 6;
                            count = count + 1;
                       else
                            count = 1; 
                       end
                   else    
                 
                       if count==1
                          MzShim(i,j) = 0.5*Gc(i,1).*te(j)*1E3;
                       elseif count==2 
                           MzShim(i,j) = 1*Gc(i,1).*te(j)*1E3 - MzShim(i,j-1);
                       elseif count==3
                           MzShim(i,j) = 1.5*Gc(i,1).*te(j)*1E3 -  MzShim(i,j-1) -  MzShim(i,j-2);
                       elseif count==4
                            MzShim(i,j) = -MzShim(i,j-1) -MzShim(i,j-2) -MzShim(i,j-3);
                       elseif count==5
                            MzShim(i,j) = 0.5*Gc(i,2).*te(j)*1E3;
                       elseif count==6 
                            MzShim(i,j) = 1*Gc(i,2).*te(j)*1E3 - MzShim(i,j-1);
                       elseif count==7
                            MzShim(i,j) = 1.5*Gc(i,2).*te(j)*1E3  -  MzShim(i,j-1) -  MzShim(i,j-2);
                       elseif count==8
                            MzShim(i,j) = -MzShim(i,j-1) -MzShim(i,j-2)- MzShim(i,j-3); 
                       end

                       if count ~= 8;
                         count = count + 1;
                       else
                         count = 1; 
                       end
                       
                   end

               end
             end
             
             
             case  'Gc_split_max'
    
    
                 for i=1:Nslices
                       count = 1; 
                       for j=start_echo+1:Nte
                           %check if a compensation gradient is 0
                           if sum(Gc(i,:)==0)  == 1
                               idx = find(Gc(i,:) ~= 0); 
                               if count==1
                                  MzShim(i,j) = 0.2*Gc(i,idx).*te(j)*1E3;
                               elseif count==2 
                                   MzShim(i,j) = 0.4*Gc(i,idx).*te(j)*1E3 - MzShim(i,j-1);
                               elseif count==3
                                   MzShim(i,j) = 0.6*Gc(i,idx).*te(j)*1E3 -  MzShim(i,j-1) -  MzShim(i,j-2);
                               elseif count==4
                                     MzShim(i,j) = 0.8*Gc(i,idx).*te(j)*1E3 -  MzShim(i,j-1) -  MzShim(i,j-2)- MzShim(i,j-3);
                               elseif count==5
                                     MzShim(i,j) = 1*Gc(i,idx).*te(j)*1E3 -  MzShim(i,j-1) -  MzShim(i,j-2)- MzShim(i,j-3) - MzShim(i,j-4);
                               elseif count==6 
                                     MzShim(i,j) = -MzShim(i,j-1) -MzShim(i,j-2)- MzShim(i,j-3) -  MzShim(i,j-4) - MzShim(i,j-5);  
                               end

                               if count ~= 6;
                                    count = count + 1
                               else
                                    count = 1; 
                               end
                           else    

                               if count==1
                                  MzShim(i,j) = 0.33*Gc(i,1).*te(j)*1E3;
                               elseif count==2 
                                   MzShim(i,j) = 0.66*Gc(i,1).*te(j)*1E3 - MzShim(i,j-1);
                               elseif count==3
                                   MzShim(i,j) = 1*Gc(i,1).*te(j)*1E3 -  MzShim(i,j-1) -  MzShim(i,j-2);
                               elseif count==4
                                    MzShim(i,j) = -MzShim(i,j-1) -MzShim(i,j-2) -MzShim(i,j-3);
                               elseif count==5
                                    MzShim(i,j) = 0.33*Gc(i,2).*te(j)*1E3;
                               elseif count==6 
                                    MzShim(i,j) = 0.66*Gc(i,2).*te(j)*1E3 - MzShim(i,j-1);
                               elseif count==7
                                    MzShim(i,j) = 1*Gc(i,2).*te(j)*1E3  -  MzShim(i,j-1) -  MzShim(i,j-2);
                               elseif count==8
                                    MzShim(i,j) = -MzShim(i,j-1) -MzShim(i,j-2)- MzShim(i,j-3); 
                               end

                               if count ~= 8;
                                 count = count + 1
                               else
                                 count = 1; 
                               end

                           end

                       end
                     end
    
    end
    
   GzShim =  cumsum(MzShim,2)./repmat(te*1E3, [Nslices, 1]); %mT/m

end

