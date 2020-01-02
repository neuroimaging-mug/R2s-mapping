function [flipAngleMap, B1Map] = B1_mapping_BS(PathName, FileName1, pulse, deltaf, alphaBS)

gamma = 2*pi*42.58e6;           %rad/s/T
deltat = 1e-6;                  %s
deltaOmega = 2*pi*deltaf;         %rad/s

info = dicominfo([PathName, FileName1]);
bitDepth = info.BitDepth;

start_idx = str2double(FileName1(end - 4));
FileName2 = [FileName1(1:end-6), num2str(1+start_idx, '%02.0f'), FileName1(end-3:end)];


phaseImagePos = double(dicomread([PathName, FileName1]));
phaseImageNeg = double(dicomread([PathName, FileName2]));

phaseImagePos = (phaseImagePos - double(2^(bitDepth-1))) / double(2^(bitDepth)) * 2*pi;
phaseImageNeg = (phaseImageNeg - double(2^(bitDepth-1))) / double(2^(bitDepth)) * 2*pi;

phaseImage = phaseImageNeg - phaseImagePos;
phaseImage(phaseImage<0) = phaseImage(phaseImage<0) + 2*pi;

Kbs = sum( pulse.^2) * gamma^2 * deltat / (deltaOmega) / (max(pulse))^2;
B1Map = sqrt(phaseImage / Kbs);
flipAngleMap = B1Map * gamma*sum(pulse)*deltat/max(pulse) *180/pi / alphaBS *100;