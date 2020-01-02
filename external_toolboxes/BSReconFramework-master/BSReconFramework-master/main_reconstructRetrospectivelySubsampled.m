%% add support functions to path
clear all
close all
clc

addpath([pwd,filesep,'supportFunctions']);
addpath([pwd,filesep,'create_mat_file']);
addpath([pwd,filesep,'nabla']);
addpath([pwd,filesep,'tgv']);
addpath([pwd,filesep,'gpu']);
setenv('LD_LIBRARY_PATH','/usr/lib/x86_64-linux-gnu')

%% define parameters
lambda = 70;                    %Regularizaiton parameters
mu = 5e-4;
implementation = 'GPU';         %'GPU': use GPU implementation
                                %'CPU': use CPU implementation
patternType = 'block';          %specify pattern type
                                %   'full':     fully sampled
                                %   'block':    block pattern in k-space center
                                %   'eliptic':  eliptical pattern in k-space center
                                %   'vdrandom': variable density pattern
                                %   'gauss':    pattern with Gaussian density function
                                
maxitTGV = 1000;                %maximum numbers of iteration                                
maxitH1  = 1000;


%% load data
FileName = 'gre_BlochSiegert_3D.dat';
PathName =  './data/';
inFile = [PathName, FileName]; 
coilSens = 'calcFromFullData';      %Coil Sensitivity Maps:
                                    %   'calcFromFullData': coil sensitivities are calculated ouf 
                                    %                       of fully sampled data
                                    %   'loadExternal':     coil sensitivities are provided externally
CoilSensPath = '';
                                 
% format of data3D: [NCol, NLine, NSlice, NPh, NCh]
[raw_data, hdr, PhaseEncDir] = readRAW_BS_3D(inFile);
raw_data = double(raw_data);
disp('Data read completed');

numSlice = hdr.Config.NImagePar;

if strcmp(PhaseEncDir,'LIN')
    dataSize = [hdr.Config.NLinMeas, hdr.Config.NImageCols, hdr.Config.NParMeas];
    imgRes = [hdr.Config.PeFOV, hdr.Config.ReadFoV, hdr.MeasYaps.sSliceArray.asSlice{1}.dThickness]./[dataSize(1:2), numSlice];
elseif strcmp(PhaseEncDir, 'COL')
    dataSize = [hdr.Config.NImageCols, hdr.Config.NLinMeas, hdr.Config.NParMeas];
    imgRes = [hdr.Config.ReadFoV, hdr.Config.PeFOV, hdr.MeasYaps.sSliceArray.asSlice{1}.dThickness]./[dataSize(1:2), numSlice];
end

if max(dataSize ~= [size(raw_data,1), size(raw_data,2), size(raw_data,3)])
    error ('Dimension missmatch')
end


%% Load Parameters
par = loadParameters();
par.impl = implementation;

%% generate pattern

switch patternType
        
    case 'full'
        pattern = ones(size(raw_data));
        
    case 'block'
        Ny = 12;
        Nz = 4;
        pattern = getRectPattern(Ny,Nz, size(raw_data));
        
    case 'eliptic'
        Ny = 12;
        Nz = 4;
        pattern = getElipticPattern(Ny,Nz, size(raw_data));
        
    case 'vdrandom'
        AF = 85; p = 14.4;
        pattern = getVdrandomPattern(size(raw_data), AF, p, numSlice);
        
    case 'gauss'
        AF = 85;
        sigma_y = 5; sigma_z = 2;
        pattern = genGaussPattern(size(raw_data), AF, [sigma_y, sigma_z], numSlice);
        
    otherwise
        error('Pattern not defined!')
        
end

par.pattern = logical(pattern(:,:,:,:,1));


%% struct par init
[dimY, dimX, dimSlice, Nph, NC] = size(raw_data);
par.NC                = NC;                       % number of coils
par.dimY              = dimY;                     % data dimensions in y
par.dimX              = dimX;                     % data dimensions in x
par.dimSlice          = dimSlice;                 % data dimensions in z

par.y                  = raw_data.*pattern;
par.PhaseEncDir        = PhaseEncDir;
scaleDim               = max(imgRes(1:2));
par.dx                 = imgRes(1)/scaleDim;
par.dy                 = imgRes(2)/scaleDim;
par.dz                 = imgRes(3)/scaleDim;

par.maxitTGV = maxitTGV;
par.maxitH1  = maxitH1;
     
par.traj='cart';            %only cartesian implemented
par.coilSens = coilSens;
par.CoilSensPath = CoilSensPath; 
if strcmp(par.coilSens, 'calcFromFullData')
    par.CoilSensData = squeeze(raw_data(:,:,:,1,:));
end

clear pattern
disp('Initialization completed');toc


%% starte reco
par.lambda = lambda;
par.mu = mu;

[B1Map, flipAngleMap]= BlochSiegReco(par);

%% remove slice oversampling
sliceDiff = (hdr.Config.NParMeas - numSlice)/2;
flipAngleMap = flipAngleMap(:,:,sliceDiff+1:end-sliceDiff);
B1Map = B1Map(:,:,sliceDiff+1:end-sliceDiff);
pattern = par.pattern(:,:,sliceDiff+1:end-sliceDiff,:,1);

%% write results into data folder
% output: 		flipAngleMap: 	relative flip angle map in % of the nominal flip angle	
%               B1Map:          B1+ map in µT
%               pattern:        used subsampling pattern

save([PathName, '/B1Map_retrospectivelySubsampled.mat'], 'flipAngleMap','B1Map', 'pattern')


