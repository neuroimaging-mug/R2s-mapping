%Script tries to implement the MERA approach for mGRE data


 
gamma = 267.51;  

te = [0:5:200]*1E-3; %s

T1exin = 0.5*1200; %ms
T1my = 10; %ms
Amy = 0.2; 
Aexin = 0.8; 
T2my = 10; %ms
T2exin = 60; %ms

%add macroscopic signal decay. 
Gz = 0.3; 
z0 = 4; 
Fn = abs(sinc(gamma*Gz*te/2)'); 
Fn = [ones(1,5) repmat([0.0, 0.2, 0.5, 0.7,0.8, 1], [1,6])]'; 

MWF_true = Amy/(Amy + Aexin)

Smy = Amy*exp(-te*1E3./T2my); 
Sexin = Aexin*exp(-te*1E3./T2exin); 
S = (Smy + Sexin).*abs(Fn)';


% Fn(:) = 1; 
fitting.Fn = Fn; 


figure; 
plot(te, S); 


data.t = te'; 
data.D = 30*S'; 





fitting.regtyp = 'me';
% % optional inputs for 'mg' fitting
% fitting.widthgauss = 5;
% fitting.numbergauss = 2;
% fitting.T0gauss = [0.02 0.07]';
% fitting.numberT0 = 5;

% for conventional regularization, the regulization weighting must be
% adjusted. This can be done manually (fitting.regadj = 'manual') or using
% one of two automatic methods ('gcv' and 'erinc') -- see literature
% references in code
fitting.regadj = 'manually';
 fitting.regweight = 0.1;

% graph the results or not. This input is irrelevant if
% analysis.interactive = 'y'

analysis.graph = 'n';

% define the range of time constants to fit. Note that for 'mg' fitting,
% this is the full range of the spectrum, but the lowest and highest mean
% component time constants (echoed to display) cover a narrower domain
fitting.rangeT = [2.5e-3 0.25];

% define the number of time constants in the spectrum
fitting.numberT = 200;

% set the non-negative least square (NNLS) code. In most cases, the supplied
% nnlsmex fuction is the fastest NNLS solver. You may need to compile this
% for your system. Or you can specify the MATLAB fuction,
% fitting.nnlscode = 'lsqnonneg'. 

% You can automatically or manually extract a finite number of component
% amplitudes and time constants. Test this out with the interactive GUI.
analysis.extract = 'auto';
analysis.numberextract = 2;




tic
[out1D,fittingout]=MERA_with_Fm(data,fitting,Fn,analysis);
toc
t_spec = out1D.T;
[val, idx] = min(abs( out1D.T - 0.025)); 
 MWF = sum(out1D.S(1:idx,:))./sum(out1D.S)
figure; plot(out1D.T, out1D.S)

% A = [23 42 37 15 52];
% M = min(A)





