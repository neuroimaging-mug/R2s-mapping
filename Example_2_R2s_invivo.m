%% Example of R2s esimtation in presence of macroscopic field variations. 
%Script performs the steps for R2s estimation as described in for single 
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


addpath('external_toolboxes'); 
addpath('GRE_modelling_toolbox'); 

%% Load parameters of a single subject


[ par_c ] = loadExampleSubjectInfoCell(  )

%% Start R2* estimation 

data = []; 
for i=1:length(par_c)
    disp(['R2s estimation: ', par_c{i}.meas_id]);
    par_tmp = par_c{i};
    [ tmp ] = prepareAndPerformR2sEstimation( par_tmp);
     data = setfield(data, par_c{i}.meas_id, tmp); 
end

%%
% noB1Slc_Gphi0z
% acq_short = {'GsPos', 'GsNeg'}; 
% sbj = fieldnames(data);
% flip_angles{1} = 'alpha_85deg'; 
% navi = {'navi_off', 'navi_on'}; 
% SE = strel('disk',3); 
% for i=1:length(sbj); 
%    for j=1:length(par.acquisitions); 
%        for k=1:length(navi); 
%            R2s{i,j,k,1} = imrotate(getfield(data, sbj{i}, par.acquisitions{j},  'results', navi{k}, flip_angles{1}, 'R2s_corr_BlochNoB1NoSlcCorr'),90);  
% 
% 
%            R2s_name{i,j,k,1} = [sbj{i}, ' ',  acq_short{j}, '' , navi{k}, ' R2s_corr_BlochNoB1NoSlcCorr ']; 
% 
%        end
%    end
% end
% 
% 



