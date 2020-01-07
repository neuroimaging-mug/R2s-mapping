%% Example of MWF esimtation in presence of macroscopic field variations. 
%Script performs the steps for MWF estimation as described in for single 
%subject: 
%
%   Soellradl M, Lesch A, Strasser J, et al. 
%   Assessment and correction of macroscopic field variations in 2D spoiled 
%   gradient-echo sequences. Magn Reson Med. 2019;00:114
%   https ://doi.org/10.1002/mrm.28139
%
% Author: Martin Soellradl
% Department of Neurology, Medical University of Graz, Graz, Austria
% email:martin.soellradl@medunigraz.at
% Website: http://www.neuroimaging.com
% Januray 2020; Last revision: 02-Januray-2020

%% 

addpath(genpath('external_toolboxes')); 
addpath(genpath('GRE_modelling_toolbox')); 


%%

[ par_c ] = loadExampleSubjectInfoCellMWF(  )


%% Start MWF estimation


data = []; 
%subjects for evaulation in paper: [1,3,4,5,7]
for i= 1:length(par_c)
    disp(['MWF estimation: ', par_c{i}.meas_id]);
    par = par_c{i} 
    par.beta_reg = 0.1;
    [ tmp ] = prepareAndPerformMWFEstimation_v2(par);
    
     data = setfield(data, par_c{i}.meas_id, tmp); 
end

close all;


%%

sbj = fieldnames(data);
flip_angles{1} = 'alpha_85deg'; 
navi = {'navi_on'}; 
SE = strel('disk',3); 
for i=1:length(sbj); 
   mask = imrotate(getfield(data, sbj{i}, par.acquisition,  'results', navi{1}, flip_angles{1}, 'mask'),90);  
   mask_er = imerode(mask, SE); 
   for j=1:length(navi); 
       MWF{i,j,1} = imrotate(getfield(data, sbj{i}, par.acquisition,  'results', navi{j}, flip_angles{1}, 'MWF_MERA_map'),90).*mask_er;  
       tmp = getfield(data, sbj{i}, par.acquisition,  'results', navi{j}, flip_angles{1}, 'MWF_MERA_noB1_noSlcCorr');
       MWF{i,j,2} = imrotate(tmp{1},90).*mask_er;  
       
      
       tmp = getfield(data, sbj{i}, par.acquisition,  'results', navi{j}, flip_angles{1}, 'MWF_MERA_wB1_wSlcCorr');
       MWF{i,j,3} = imrotate(tmp{1},90).*mask_er;  
       
       MWF_name{i,j,1} = [sbj{i}, ' ', navi{j}, ' MWF_MERA_map ']; 
       MWF_name{i,j,2} = [sbj{i}, ' ', navi{j}, ' MWF_MERA_noB1_noSlcCorr ']; 
       MWF_name{i,j,3} = [sbj{i}, ' ', navi{j}, ' MWF_MERA_wB1_wSlcCorr ']; 
       
       mag{i,j} = imrotate(getfield(data, sbj{i}, par.acquisition,  'results', navi{j}, flip_angles{1}, 'mag'),90);  
       mag_name{i,j} = [sbj{i}, ' ', navi{j}, ' mag ']; 
       
       Gz{i,j} = imrotate(getfield(data, sbj{i}, par.acquisition,  'results', navi{j}, flip_angles{1},'Gz'), 90); 
       Gz_name{i,j} = [sbj{i}, ' ', navi{j}, ' Gz '];
   end
end

%% Show results with imagine 

vx_size = [1.1,1.1,4]; 
plotCellImagine(squeeze(MWF), squeeze(MWF_name), vx_size); 

