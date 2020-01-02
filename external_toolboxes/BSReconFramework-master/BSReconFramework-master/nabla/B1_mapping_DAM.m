function [flipAngleMap] = B1_mapping_DAM(PathName, FileName)%, pulse)

% gamma = 2*pi*42.58e6;           %rad/s/T
% deltat = 1e-6;                  %s

inFile = fullfile(PathName,FileName);
info = dicominfo(inFile);
alpha = info.FlipAngle;
SI60 = double(dicomread(inFile));

start_idx = str2double(FileName(11:12));
filename2 = [FileName(1:10), num2str(1+start_idx, '%02.0f'), FileName(13:end)];
SI120 = double(dicomread([PathName, filename2]));

alp60=acos(SI120./(2*SI60))*180/pi;
alp60(isnan(alp60))=0;alp60(isinf(abs(alp60)))=0;alp60(angle(alp60)~=0)=0;%alp60(alp60==0)=alpha;

flipAngleMap = alp60/alpha*100;
%alp60 = alp60*pi/180;

%B1Map = alp60 / gamma / (sum(pulse)*deltat) * max(pulse);
