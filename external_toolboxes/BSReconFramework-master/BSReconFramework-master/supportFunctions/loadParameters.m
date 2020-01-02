function par = loadParameters()

file1 = fopen('BS_pulse.txt');
if file1 == -1
    error(['File ', 'BS_pulse.txt',' does not exist!']);
else
    fclose(file1);
end

data = importfile('BS_pulse.txt');
par.pulse = data(data>0);

fileID = fopen('parameters.txt');

if fileID == -1
    error(['File ', 'parameters.txt',' does not exist!']);
end

C = textscan(fileID, '%s');
par.Tpulse = str2double(C{1}{2});
par.deltaf = str2double(C{1}{5})*1000;
par.deltaOmega = 2*pi*par.deltaf;
par.alphaBS = str2double(C{1}{8});
par.gamma = 2*pi*42.58e6;           %rad/s/T
par.deltat = 1e-6;                  %s
par.Kbs = sum( par.pulse.^2) * par.gamma^2 * par.deltat / (par.deltaOmega) / (max(par.pulse))^2;
fclose(fileID);