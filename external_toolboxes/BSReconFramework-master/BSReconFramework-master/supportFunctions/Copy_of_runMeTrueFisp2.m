%% add support functions to path
clear all
close all
clc
addpath([pwd,filesep,'supportFunctions']);
addpath([pwd,filesep,'matlabCgDescent']);
addpath([pwd,filesep,'exp-MARTINI-T1']);
addpath([pwd,filesep,'create_mat_file']);
addpath([pwd,filesep,'slice_functions']);
addpath([pwd,filesep,'referenceFunctions']);
 
%% lade Phantom daten
threshold = input('Please type in the threshold for generating the mask:\n');
radius = input('Please type in radius for mask closing operation:\n');
[FileName_raw,PathName_raw] = uigetfile('*.dat','Select the raw-data');
inFile_raw = fullfile(PathName_raw,FileName_raw);

[FileName_dicom,PathName_dicom] = uigetfile('*.DCM','Select the first file of the TRUFI sequence',PathName_raw);
inFile_dicom = fullfile(PathName_dicom,FileName_dicom);

[FileName_DAM,PathName_DAM] = uigetfile('*.DCM','Select the first file of the DAM sequence for B1 correction',PathName_dicom);
inFile_DAM = fullfile(PathName_DAM,FileName_DAM);

[data, N_seg] = readRAW_multichannel(inFile_raw);

info = dicominfo(inFile_dicom);
N_ramp = 20;
FA = info.FlipAngle;
N_phase = info.CardiacNumberOfImages;
TR_pulse = info.RepetitionTime/1000/N_seg;
bitdepth = double(info.BitDepth);
TT = info.TriggerTime/1000;
TR = info.RepetitionTime/1000;
TI = info.TriggerTime/1000 - N_ramp*TR_pulse;
Tmax = TR_pulse * N_seg * N_phase;
Nmax = N_seg * N_phase;
alpha0 = 180*pi/180;
precond_mod = 2;
% threshold = 0.11;
undersampling_mode = 3;
maxit = 200;

if inFile_DAM
    corr_factor = B1_corr(PathName_DAM, FileName_DAM);
    corr_factor = imresize(corr_factor, [info.Height, info.Width], 'bicubic');
else
    corr_factor = ones(info.Height, info.Width);
end

% load raw_Data.mat
% TR = 36.64/1000;    % in s
% TT = 115/1000;      % in s
% FA = 40;            % in deg

% TI = 0.0234;
% N_seg = 8;
% N_phase = 89;
% TR_pulse = 4.58/1000;


%% Undersampling Mode
switch undersampling_mode
    case 1
        % simuliere standard interleaved undersampling
        AF        = 4;  % simulated accerelation factor
        ACL       = 0;  % autocalibration lines
        symmetric = false;
        
        uData      = interleavedPattern(data,AF,ACL,symmetric);
        
    case 2
        % ... alternativ: simuliere blocked undersampling
        AF       = 3;              % simulated accerelation factor
        ordering = 'linear';
        startAt  = 'centerblock';
        
        uData = blockedPattern(data,'AF',AF,'ordering',ordering,'startat',startAt);
        
    case 3
        % ...oder nimm halt erst mal vollgesamplete Daten :)
        uData = double(data(:, :, : , :));
end

%% struct par init
clear data
[dimY, dimX, NPh, NC] = size(uData);
par.TR                = TR;                       % repetition time
par.NPh               = NPh;                      % number of phases
par.NC                = NC;                       % number of coils
par.dimY              = dimY;                     % data dimensions
par.dimX              = dimX;                     % data dimensions
par.dims              = [dimY, dimX];             % data dimensions
par.fa                = (FA / 180) * pi;          % flip angle in rad
par.TI                = TI;
par.Nseg              = N_seg;
par.TR_pulse          = TR_pulse;
par.Nmax              = Nmax;
par.Tmax              = Tmax;
par.alpha0            = alpha0;
par.precond_mod       = precond_mod;
par.Nramp             = N_ramp;
par.corr_factor       = corr_factor;

% 'standardize the data'
par.yscale = 1 / sqrt(scalar(uData(:))) * NC;
par.y      = uData * par.yscale;

% CG parameter
par.maxit               = maxit;      % maximum number of iterations
par.tol                 = 1e-10;    % CG tolerance
par.cgRestarts          = 1;        % number of CG restarts (after every 
                                    % restart, the scaling parameters are updated)
                                    
% coil sensitivity estimation
if NC==1
    par.coilEstimate        = 'ones';      % can be either 'nlinv' or 'ones'
                                           % (see below)
else 
    par.coilEstimate        = 'nlinv';
end

%% starte reco (fuer groessere Daten ist das bisher leider ziemlich zeitaufwaendig)
clear uData info corr_factor;
[result,  guess, par] = trueFispReco(par, FA, TT, threshold, radius);

fig1 = figure;
imshow(result(:,:,1)*1000,[0 4000])
fig2 = figure;
imshow(result(:,:,2)*1000,[0 400])
fig3 = figure;
imshow(result(:,:,3),[0,1])
fig4 = figure;
imshow(par.mask,[0 1])

outdir = datestr(now, 'yyyy-mm-dd HH_MM_SS');
mkdir('./output', outdir);

saveas(fig1, ['./output/', outdir, '/T1_map.fig']);
saveas(fig2, ['./output/', outdir, '/T2_map.fig']);
saveas(fig3, ['./output/', outdir, '/M0_map.fig']);
saveas(fig4, ['./output/', outdir, '/mask.fig']);

save(['./output/', outdir, '/output_data.mat'], 'result', 'guess', 'par');

