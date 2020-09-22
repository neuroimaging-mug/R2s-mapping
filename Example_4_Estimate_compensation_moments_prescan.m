%% Example for estimating the slice-specific field gradient values from the 
%prescan data. 

%Script performs the steps for R2s estimation as described in:
%   Soellradl, M, Strasser, J, Lesch, A, Stollberger, R, Ropele, S, Langkammer, C. 
%   Adaptive slice-specific z-shimming for 2D spoiled gradient-echo sequences. 
%   Magn Reson Med. 2020; 00: 1-14. https://doi.org/10.1002/mrm.28468 

% Author: Martin Soellradl
% Department of Neurology, Medical University of Graz, Graz, Austria
% email:martin.soellradl@medunigraz.at
% Website: http://www.neuroimaging.com
% Januray 2020; Last revision: 02-Januray-20

path_pwd = pwd; 

%% Before running the example please make sure that you have dowloaded the 
%provided input data from zenodo: 
 
%https://doi.org/10.5281/zenodo.4044429
 
% Please download data and unzip it in the main folder. The folder contains
% also the results. 

if ~exist([pwd, '/data_input'], 'dir'); 
   error('Please download input data from https://doi.org/10.5281/zenodo.4044429 and extract it into the current directory'); 
end



%% 

addpath(genpath('external_toolboxes')); 
addpath(genpath('GRE_modelling_toolbox')); 


%%

% gyromagnetic ratio 42.57*2*pi in rad*1E6/T
gamma = 267.51;           
    
meas_id = 'subject_1';

src_nii =  [path_pwd, '/data_input/Example_4_zshimTable_prescan/', meas_id, '/nifti/']
path_results = [path_pwd, '/results/R2s_zshim/', meas_id, '/prescan/']; 

if ~exist(path_results, 'dir')
    mkdir(path_results); 
end


%% Load data 



tmp = load_untouch_nii([src_nii, 'zShim_prescan_phase', '.nii.gz']); 
phase_siemens = double(tmp.img); 

%convert to 2pi 
ma = 4094;
mi =  -4096;
dif = ma-mi;
phase  = (phase_siemens-mi)/(dif)*2*pi-pi;


tmp = load_untouch_nii([src_nii, 'zShim_prescan_mag', '.nii.gz']); 
mag = double(tmp.img); 

%store header of the nifit files 
tmp.hdr.dime.datatype = 16; 
tmp.hdr.dime.bitpix = 32;   
nii_template_4D = tmp;
tmp.hdr.dime.dim([1 5])=[4 1];
nii_template = tmp; 


%slice thickness
z0 = 3; %mmm

%Echo times
te = [2.69, 4.81, 6.93]; %ms


%% Estimate the field gradient map from the pre-scan data

opts = []; 
opts.nii_template = nii_template;
opts.nii_template_4D = nii_template_4D;
opts.B0_method  = 'phase_diff_B0';
opts.bbet = 1; 
opts.bprelude = 1;
file_id = 'ref_scan'; 

[ res] = estimateGradientMaps( mag, phase, te, file_id, path_results, opts );


%% Post processing steps of the field gradient map

%erode mask 
mask = res.mask;
SE = strel('disk',5);
mask_er = imerode(mask, SE); 



%Low pass filter Gz
Gz =  imgaussfilt3(1E3*res.Gz.*mask_er./(gamma*z0),1); 

Gz_sub = zeros(size(Gz)); 
th = [-500:10:500]; 
for i=1:size(Gz,3); 
    for k=1:length(th)-1; 
       th_vec = [th(k), th(k+1)]; 
       [ mask_Gz ] = thresholdImage(Gz(:,:,i),th_vec );
      
       th_mean(k) = mean(th_vec);
       Gz_sub(:,:,i) = Gz_sub(:,:,i) + mask_Gz*mean(th_vec); 

    end
end
Gz_sub = Gz_sub.*mask_er; 

Gz_c{1} = Gz_sub; 
mask_c{1} = mask_er;
name_c{1}= 'Gz'; 


%% Calculate compenstation moment (min/max values in each slice)

Gz_sep = []; 
for i=1:size(Gz,3); 
   Gz_tmp = Gz_sub(:,:,i); 
   Gz_vec =Gz_tmp((mask_er(:,:,i) > 0));
   
   Gz_low = min(Gz_vec(Gz_vec < 0)); 
   Gz_high= max(Gz_vec(Gz_vec > 0));  

   if isempty(Gz_low)
        Gz_low = 0; 
   end
   
   if isempty(Gz_high)
        Gz_high = 0; 
   end
   Gz_sep(i,:) =[Gz_low, Gz_high]; 
    
    
end


Gz_comp = -1E-3*Gz_sep; %mT/m
Gz_comp(abs(Gz_comp) > 25000) = 0; 


%% Write Moments into table 

filename = [path_results, 'zShimTable_', meas_id, '.txt']; 
fid = fopen(filename,'wt');
 for ii = 1:size(Gz_comp,1)
    fprintf(fid,'%d %d\r\n',Gz_comp(ii,:));
end
fclose(fid)

