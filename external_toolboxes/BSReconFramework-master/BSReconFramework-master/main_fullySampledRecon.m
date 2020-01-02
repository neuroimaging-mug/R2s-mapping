%% add support functions to path
clear all
close all
clc

addpath([pwd,filesep,'supportFunctions']);
addpath([pwd,filesep,'create_mat_file']);
addpath([pwd,filesep,'nabla']);
addpath([pwd,filesep,'tgv']);

%% load data
FileName = 'meas_MID00200_FID30180_gre_BlochSiegert_3D.dat';
PathName =  './data/';
inFile = [PathName, FileName];

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

%% load parameters
par = loadParameters();

%% start coil combination 

disp('Starting Walsh for positive full Dataset')
dataPos = squeeze(raw_data(:,:,:,1,:));
imgPos = walsh_sens_3d_slice(ifft3c(dataPos));
toc

disp('Starting Walsh for negative full Dataset')
dataNeg = squeeze(raw_data(:,:,:,2,:));
imgNeg = walsh_sens_3d_slice(ifft3c(dataNeg));
toc

%% Calculate B1+ field
disp('Calculating Bloch Siegert Map for full Dataset')
phaseImage = angle(imgNeg) - angle(imgPos);
phaseImage(phaseImage<0) = phaseImage(phaseImage<0) + 2*pi;
B1Map_full = sqrt(phaseImage / par.Kbs);
flipAngleMap_full = B1Map_full * par.gamma*sum(par.pulse)*par.deltat/max(par.pulse) *180/pi / par.alphaBS *100;

sliceDiff = (hdr.Config.NParMeas - numSlice)/2;
flipAngleMap_full = flipAngleMap_full(:,:,sliceDiff+1:end-sliceDiff);
B1Map_full = B1Map_full(:,:,sliceDiff+1:end-sliceDiff);

save('./data/fullySampledReference.mat', 'B1Map_full', 'flipAngleMap_full')

