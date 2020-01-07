%% Example of R2s esimtation in presence of macroscopic field variations for a 
%single subject. 
%Script performs the steps for R2s estimation as described in:
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

%% Before running the example please make sure that you have dowloaded the 
%provided input data from zenodo: 
 
%https://doi.org/10.5281/zenodo.3600319
 
% Please download data and unzip it in the main folder. Also, the results
% are available, which allows to just load the results. 

if ~exist([pwd, '/data_input'], 'dir'); 
   error('Please download input data from https://doi.org/10.5281/zenodo.3600319 and extract it into the current directory'); 
end


%% Set paths

addpath(genpath('external_toolboxes')); 
addpath(genpath('GRE_modelling_toolbox')); 

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



