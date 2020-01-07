function [  ] = plot4DCellImagine( data, names,  voxel_size, opt, Pos,  te)
%PLOTCELLIMAGINE Plots an cell array that contains 4D volumes. The name for
%each volume is stored in the names cell array. For each volume either one
%slice position in z (3 dimension) or in time domain (4 dimension) is
%plotted with imagine. 
%   @param      data     cell array containting 4D 
%   @param      name     cell array with names of for each 4D volume  
%   @param      opt      String opt for plotting either 3D for specific
%                        echo numer of as a function of echo time for one
%                        slice position in z. 
%                        opt = 'SlicePos' or opt = 'EchoNumber'
%   @param      te       echo times in ms (optional) 


    var_name = 'data'; 
    all_img = []; 
    for i=1:1:size(data,1); 
       for j=1:size(data,2); 
          for  k=1:1:size(data,3); 
             
              

            switch opt
                case 'SlicePos' 
                      img_call = ['squeeze(', var_name, '{', num2str(i), ',', num2str(j), ',', num2str(k), '}(:,:,', num2str(Pos), ',:))']; 
                case 'EchoNumber'
                      img_call = ['squeeze(', var_name, '{', num2str(i), ',', num2str(j), ',', num2str(k), '}(:,:,:,', num2str(Pos) '))'];
                otherwise
                    error('Invalid input for opt. Must be string containing either the option ''SlicePos'' or ''EchoNumber''')
            end
           
             par_vx = ['''VoxelSize''', ', [', num2str(voxel_size), ']']; 
             name_img = ['''Name''', ',''', names{i,j,k}, '''' ];

             sgl_img = [img_call, ',' par_vx, ' , ' name_img, ',']; 

             all_img = [all_img, sgl_img]; 
          end
       end
    end
    
    if nargin == 5
        all_img(end) = []
        imagine_call = ['imagine(', all_img, ')']; 
    else
        imagine_call = ['imagine(', all_img, '''te''', ', te', ')']; 
    end
    eval(imagine_call); 
end

