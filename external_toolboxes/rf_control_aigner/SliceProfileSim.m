%Script that simulates different pulse B1 pulse shape. For the forward
%simluation the CN Bloch equations are used. 


%First define the timings of the pulse and slice selection gradient 
Tpulse = 2;                  %Pulse duration in ms 
Tramp = 0.064;                %Ramp up/down time of Gz in ms
Tsim = 4*Tramp + Tpulse + Tpulse/2; 

d.T     = Tsim;              % optimization time in ms 
d.dt = 2E-3;               % temporal grid size in ms 
d.Nt = Tsim/d.dt;            % total number of temporal points
d.tdis = 0:d.dt:Tsim-d.dt;   % temporal running variable

%calculate B1 (gauss)
tB1max = Tpulse/2 + Tramp;   %Maximum of B1 
Amp = 0.8367/2;                     %weird unit thing [Amp] = 0.1 µT
sigma = 0.280;                  %ms 
B1 = Amp*exp(-(d.tdis - tB1max).^2./(2*sigma.^2)); 
% B1 = Amp*sin(2*pi*(d.tdis + tB1max)./Tpulse); 

pulse_sample = round(Tpulse/d.dt); 
% t_start = 0.1; %ms
% Amp = 1.156;
% B1 = zeros(1, d.Nt)
% B1(t_start/d.dt: t_start/d.dt + pulse_sample) = Amp; 

%estimate flip angle 

d.gamma = 267.51;       % gyromagnetic ratio 42.57*2*pi in MHz/T
alpha = rad2deg(d.gamma.*trapz(d.tdis*1E-3, B1*10))


%% Calculate slice selective gradient 


Gz_amp = 6.58;                 %plateau amplitude of Gz in mT/m
Gz = zeros(1,d.Nt);

%ramp down 
ramp_end = round(Tramp/d.dt);
ramp_sample = 1:ramp_end; 

k = -Gz_amp./(ramp_sample(end) - ramp_sample(1));
doff = -Gz_amp - k*ramp_sample(end);
Gz(ramp_sample) =k.*ramp_sample + doff;

%plateu 
plateu_end = ramp_end + round(Tpulse/d.dt);
Gz(ramp_end:plateu_end) = -Gz_amp; 

%ramp up till rephaser 
%ramp_sample = (plateu_end):plateu_end  + 2*ramp_end -1; 
ramp_sample = (plateu_end):plateu_end  + ramp_end -1; 
k = 2*Gz_amp./(ramp_sample(end) - ramp_sample(1));
doff = Gz_amp - k*ramp_sample(end);
Gz(ramp_sample) =k.*ramp_sample + doff;


% %plateu rephaser
% plateu_sample = ramp_sample(end) : ramp_sample(end)+(Tpulse/(2*d.dt)) +1;
% Gz(plateu_sample) = Gz_amp; 
% 
% %ramp down 
% ramp_sample = plateu_sample(end): plateu_sample(end)+round(Tramp/d.dt);
% k = -Gz_amp./(ramp_sample(end) - ramp_sample(1));
% doff = Gz_amp - k*ramp_sample(1);
% Gz(ramp_sample) =k.*ramp_sample + doff;

figure; 
subplot(2,1,1); 
plot(d.tdis, B1*10);
xlabel('time in ms'); 
ylabel('B1 in µT'); 
title('RF Pulse'); 
subplot(2,1,2); 
plot(d.tdis, Gz); ylim([-1.4*Gz_amp, 1.4*Gz_amp]); 
xlabel('time in ms'); 
ylabel('Gz in mT/m'); 
title('Slice selection gradient'); 


%% set parameters of Christoph's script 



% time discretization
load('z_grad_thk2_dt5.mat');     % slice selective gradient shape


% space discretization
d.a    = 0.25;                     % domain border in m
d.z    = 0.0025;                  % half slice thickness in m
d.Nx   = 2501;                    % total number of spatial points
d.xdis = linspace(-d.a,d.a,d.Nx); % spatial running variable
d.dx   = d.xdis(2)-d.xdis(1);     % spatial grid size

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
d.v  =  zeros(d.Nt,1); % fixed with zeros
d.w  =    Gz;          % fixed with external shape
d.alpha = 1e-4 ;        % control costs for u (SAR)
d.Nu    = d.Nt;                   % number of temporal control points

u = B1; 

% initial magnetization
d.M0    = d.M0c*repmat([0;0;1],1,d.Nx); 


M = cn_bloch(d,d.M0,u,d.v,d.w ); 

% figure; 
% Mxy = sqrt(M(1, :,end).^2 + i*M(2,:,end).^2); 
% plot(d.xdis ,Mxy);

Mxy = sqrt(M(2,:,end).^2+M(1,:,end).^2); 

[rho, Mxy] = cart2pol(M(2,:,end),M(1,:,end));


% slc_th = fwhm(d.xdis*1E3, Mxy); 

% plot control variable together with slice selective gradient
figure; 
    plot(d.tdis(1:end),u*1000*d.B1c,'r');
    hold on
    grid on
    plot(d.tdis(1:end),d.v,'Color',[0 0.5 0]);
    plot(d.tdis(1:end),d.w,'k');
    hold off
    legend('B_{1,x}','B_{1,y}','G_z', 'Location', 'NorthEastOutside');
    xlabel('time in ms');
    ylabel('B_1 in \mu T');
    axis([0, d.T, -20, 20]);

figure; 
subplot(2,1,1); % plot magnetization after excitation
    plot(d.xdis*1E3,M(1,:,end),'Color',[0 0.5 0]);
    hold on
    grid on
    plot(d.xdis*1E3,M(2,:,end),'b','LineWidth',1.5);
    plot(d.xdis*1E3,M(3,:,end),'r','LineWidth',1.5);
    hold off
    legend('M_x(T)','M_y(T)','M_z(T)');
    axis([-d.a*1E3/20, d.a*1E3/20, -0.1, 1.1]);
    xlabel('distance in mm');
    ylabel('normalized magnetization');
    
subplot(2,1,2); % plot transverse magnetization after excitation (zoom)
    plot(d.xdis*1E3,Mxy ,'b');
    grid on
    legend('M_{xy}(T)');
    axis([-d.a*1E3/20, d.a*1E3/20, -0.1, 1.1]);
    xlabel('distance in mm');
    ylabel('normalized magnetization');

% dim = [.2 .5 .2 .3];
% str = ['FWHM = ', num2str(slc_th,3), 'mm'];
% annotation('textbox',dim,'String',str,'FitBoxToText','on');



%% Solve with Fourier approximation 















