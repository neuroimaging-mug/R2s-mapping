
clear alL; 

%% set parameters
% space discretization
d.a    = 0.5;                     % domain border in m
d.z    = 0.0025;                  % half slice thickness in m
d.Nx   = 5001;                    % total number of spatial points
d.xdis = linspace(-d.a,d.a,d.Nx); % spatial running variable
d.dx   = d.xdis(2)-d.xdis(1);     % spatial grid size

% time discretization
load('z_grad_thk2_dt5.mat');     % slice selective gradient shape
d.T     = 3.480;                 % optimization time in ms 
d.Nt    = size(z_grad,1)+1;      % total number of temporal points
d.tdis  = linspace(0,d.T,d.Nt);  % temporal running variable
d.dt    = d.tdis(2) - d.tdis(1); % temporal grid size
d.Nu    = 512;                   % number of temporal control points

% model parameters
d.gamma = 267.51;       % gyromagnetic ratio 42.57*2*pi in MHz/T
%d.T1    = 102;          % longitudinal relaxation time in ms
d.T1    = 1000;          % longitudinal relaxation time in ms
d.T2    = 81;           % transversal relaxation time in ms
d.B0    = 3000;         % static magnetic field strength in mT
d.M0c   = 1;            % normalized equilibrium magnetization
d.B1c   = 1e-2;         % weighting for the RF amplitude [u*1e3*d.B1c] = muT
d.G3    = 1;            % weighting for the z-Gradient in mT
d.relax = 0;            % 0=without relaxation, 1=with relaxation
d.u0 = zeros(d.Nu,1);   % RF initial guess
d.v  = zeros(d.Nt-1,1); % fixed with zeros
d.w  = z_grad;          % fixed with external shape
d.alpha = 1e-4 ;        % control costs for u (SAR)

% initial magnetization
d.M0    = d.M0c*repmat([0;0;1],1,d.Nx); 


%load previous optimization result 
load('u'); 
u_new = 1.5*gaussmf(-d.Nu/2:d.Nu/2-1, [50,0])'; 

% zero padding of pulse to readout time
u = [u; zeros(d.Nt-1-d.Nu,1)];
u_new = [u_new; zeros(d.Nt-1-d.Nu,1)];
v = zeros(size(u)); 


%MS adaptions 
% u is equal to B1x shape [Nt]
% v is equal to B1y shape [Nt]
% w is the gradient [Nt]
% M0 start magnetization 

figure; plot(d.tdis(1:end-1), u);
figure; subplot(2,1,1); 
plot(d.tdis(1:end-1), u_new);
subplot(2,1,2);  plot(d.w ); 

%create a gauss


M = cn_bloch(d,d.M0,u_new,v,d.w ); 

figure; 
plot(d.xdis , M(3, :,end)); 

 plot_results(u,d)  

figure; 
plot(d.xdis , M(3, :,end)); 

figure; 
Mxy = sqrt(M(1, :,end).^2 + i*M(2,:,end).^2); 
plot(d.xdis ,Mxy); 


figure; 
plot(d.xdis , phase); 

figure; 
plot(d.xdis , M(1, :,end) ); 




