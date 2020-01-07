function [  ] = plotCellImagine( data, names,  voxel_size, te)
%PLOTCELLIMAGINE Plots an cell array that contains 3D volumes. The name for
%each volume is stored in the names cell array. 
%   Detailed explanation goes here
   
    var_name = 'data'; 
    all_img = []; 
    for i=1:1:size(data,1); 
       for j=1:size(data,2); 
          for  k=1:1:size(data,3); 
             img_call = [var_name, '{', num2str(i), ',', num2str(j), ',', num2str(k), '}']; 
             par_vx = ['''VoxelSize''', ', [', num2str(voxel_size), ']']; 
             name_img = ['''Name''', ',''', names{i,j,k}, '''' ];

             sgl_img = [img_call, ',' par_vx, ' , ' name_img, ',']; 

             all_img = [all_img, sgl_img]; 
          end
       end
    end
    if nargin == 3
        all_img(end) = []
        imagine_call = ['imagine(', all_img, ')']; 
    else
        imagine_call = ['imagine(', all_img, '''te''', ', te', ')']; 
    end
    eval(imagine_call); 
end

