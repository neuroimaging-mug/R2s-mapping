%% Example of R2s esimtation from multi-echo gradient-echo (mGRE) with z-shim 
%compensation moments. The in-vivo data set example includes data from
%standard mGRE data, from data with constant +/- linear increasing
%compensation moments, and the proposed slice-specfic z-shim pattern. 

%Script performs the steps for R2s estimation as described in:
%   Soellradl, M, Strasser, J, Lesch, A, Stollberger, R, Ropele, S, Langkammer, C. 
%   Adaptive slice-specific z-shimming for 2D spoiled gradient-echo sequences. 
%   Magn Reson Med. 2020; 00: 1-14. https://doi.org/10.1002/mrm.28468 

% Author: Martin Soellradl
% Department of Neurology, Medical University of Graz, Graz, Austria
% email:martin.soellradl@medunigraz.at
% Website: http://www.neuroimaging.com
% Januray 2020; Last revision: 02-Januray-2020


%% Before running the example please make sure that you have dowloaded the 
%provided input data from zenodo: 
 
% https://doi.org/10.5281/zenodo.4044429
 
% Please download data and unzip it in the main folder. The folder contains
% also the results. 

if ~exist([pwd, '/data_input'], 'dir'); 
   error('Please download input data from https://doi.org/10.5281/zenodo.4044429 and extract it into the current directory'); 
end



%% 

addpath(genpath('external_toolboxes')); 
addpath(genpath('GRE_modelling_toolbox')); 


%%

gamma = 267.51;  
[par_c ] = loadExampleSubjectInfoAdaptivezShimming(  )



data = []; 
for i=1:length(par_c);   
    par_zShim = par_c{i}; 
	tmp = prepareAndPerformR2sEstimationWithzShim_v1_git( par_zShim );
    subjects{i} = par_c{i}.meas_id; 
    data = setfield(data, subjects{i}, tmp); 
end


%% Show results 

te = data.subject_1.zShim_Gc_220.results.navi_on.alpha_60deg.dcm_header.EchoTime; 
te(end) = []; 
acq_short = {'zShim off',  'Gc 230µT/m', 'Gc split'}; 
sbj = fieldnames(data);
navi = {'navi_off', 'navi_on'}; 
SE = strel('disk',3); 

acquisitions{1} = 'zShim_off';
acquisitions{2} = 'zShim_Gc_220';
acquisitions{3} = 'zShim_Gc_split';

z0 = data.subject_1.zShim_Gc_220.results.navi_on.alpha_60deg.dcm_header.SliceThickness;
for i=1:length(subjects); 
    for j=1:length(acquisitions); 
       flip_angles = fieldnames(getfield(data, subjects{i}, acquisitions{j}, 'results', navi{2})); 
       for k=1:length(flip_angles); 
           mask = imrotate(getfield(data, subjects{i},acquisitions{j},  'results', navi{2}, flip_angles{k}, 'mask'),90); 

            R2s{i,j,k,1} = 1E3*imrotate(getfield(data, subjects{i},acquisitions{j},  'results', navi{2}, flip_angles{k}, 'R2s_corr_zShimNoB1NoSlcCorr'),90).*mask;  
            R2s{i,j,k,2} = 1E3*imrotate(getfield(data, subjects{i},acquisitions{j},  'results', navi{2}, flip_angles{k}, 'R2s_corr_zShimWithB1NoSlcCorr'),90).*mask; 
            R2s{i,j,k,3} = 1E3*imrotate(getfield(data, subjects{i},acquisitions{j},  'results', navi{2}, flip_angles{k}, 'R2s_corr_zShimWithB1WithSlcCorr'),90).*mask; 
            
            
            mag{i,j,k} = imrotate(getfield(data, subjects{i},acquisitions{j},'results', navi{2}, flip_angles{k}, 'mag'),90).*repmat(mask, [1,1,1,length(te)]);
            mag_tmp= imrotate(getfield(data, subjects{i},acquisitions{j},'results', navi{2}, flip_angles{k}, 'mag'),90);  

            
            R2s_name{i,j,k,1} = [acq_short{j}, ' R2s_corr_BlochNoB1NoSlcCorr ']; 
            R2s_name{i,j,k,2} = [acq_short{j}, ' R2s_corr_zShimWithB1NoSlcCorr ']; 
            R2s_name{i,j,k,3} = [acq_short{j}, ' R2s_corr_zShimWithB1NoSlcCorr ']; 

            R2s_xz{i,j,k,1} = permute(R2s{i,j,k,1}, [1,3,2]);
            R2s_xz{i,j,k,2} = permute(R2s{i,j,k,2}, [1,3,2]);
            
            mask_c{i,j,k,1} = mask;
            
            Gz_tmp = getfield(data, subjects{i},acquisitions{j},'results', navi{2}, flip_angles{k}, 'Gz');
            Gz{i,j,k} = imrotate(1E3*Gz_tmp./(gamma*z0),90); 
       end

    end
end

%% Plot results with Imagine 


voxel_size = [1,1,3]; 
plotCellImagine( squeeze(R2s(1,:,1,1)), squeeze(R2s_name(1,:,1,1)),  voxel_size)

%
voxel_size = [1,1,3]; 
opt= 'SlicePos';
te = data.subject_1.zShim_Gc_220.results.navi_on.alpha_60deg.dcm_header.EchoTime; 
te(end) = []; 
Pos = 12; 
plot4DCellImagine(mag(1,[1:3]), R2s_name(1,[1:3]), voxel_size,  opt,Pos,  te);

