
function [Mcplx, FWHM, d] = simSliceProfileBlochEqn_v2( pulse, Gz_amp, Gmac)
%SIMSLICEPROFILE Sovles Bloch equations for the selected RF-pulse 
%equation with the Bloch solver provied by C Aigner, JMRI 2015
%   @param      pulse   Struct with pulse parameters
%   @param      Gz_amp  Plateu amplitude of slice selection Gradient [mT/m]
%   @param      Gmac    Macroscopic field gradient in z direction [mT/m] 
%                       (optional)
%
% v2: Polarity issue resolved.  Polarity is changed. Positive slice direction means that the it is
% positve in the gradient coordinate system, which is given by the
% cross-prodcut of the phase encdoing (e.g. R>>L) and the readout
% direction. Assuming that patient is head first the postive slice
% direciton is opposed to the device coordiante system where z points out
% of the scanner. 

%   Author: Martin Soellradl
%   Department of Neurology, Medical University of Graz, Graz, Austria
%   email:martin.soellradl@medunigraz.at
%   Website: http://www.neuroimaging.com
%   Januray 2020; Last revision: 02-Januray-2020


    if nargin < 3
       Gmac = 0;  
    end
    
    %First define the timings of the pulse and slice selection gradient 
    Tpulse = pulse.Tpulse;                  %pulse duration in ms 
    Tramp = 0.064;                          %ramp up/down time of Gz in ms
    d.Tramp = Tramp;
    Tsim = 4*Tramp + Tpulse + Tpulse/2;     %total Simulation time 

    d.T  = Tsim;                            %optimization time in ms 
    d.dt = 2E-3;                            %temporal grid size in ms 
    d.Nt = Tsim/d.dt;                       %total number of temporal points
    d.tdis = 0:d.dt:Tsim-d.dt;              %temporal running variable
    
    d.gamma = 267.51;                       % gyromagnetic ratio
                                            % 42.57*2*pi in MHz/T

    
    %% Calculate the B1 pulse envelope
    
    [ B1 ] = getB1PulseEnvelope( pulse, d.tdis, Tramp);

    
    %% Calculate slice selection gradient

    Gz = zeros(1,d.Nt);
    
    %Ramp down
    ns = round(Tramp/d.dt);
    start_idx = 1; 
    ramp_samples = start_idx:ns; 
    k = -Gz_amp./(ramp_samples(end) - ramp_samples(1));
    doff = -Gz_amp - k*ramp_samples(end);
    Gz(ramp_samples) =k.*ramp_samples + doff;

    %Plateu of Gz 
    ns = round(Tpulse/d.dt);
    start_idx = round(Tramp/d.dt);
    plateu_samples = start_idx:ns+start_idx; 
    Gz(plateu_samples) = -Gz_amp; 
    
    %Ramp up
    start_idx = round((Tramp + Tpulse)/d.dt);
    ns = 2*Tramp/d.dt; 
    ramp_samples = start_idx:start_idx + ns; 
    k = 2*Gz_amp./(ramp_samples(end) - ramp_samples(1));
    doff = Gz_amp - k*ramp_samples(end);
    Gz(ramp_samples) =k.*ramp_samples + doff;

    %Rephaser Gz 
    start_idx = round((Tramp + Tpulse + 2*Tramp)/d.dt);

    Tplateu = Tpulse/2 - Tramp/2; %Bernstein, p. 76. Rephasing moment AR
    ns = round(Tplateu/d.dt);
    plateu_samples = start_idx:start_idx+ns; 
    Gz(plateu_samples) = Gz_amp; 
    
    %Ramp down 
    start_idx = round((Tramp + Tpulse + 2*Tramp + Tplateu)/d.dt);
    ns = round(Tramp/d.dt); 
    ramp_samples = start_idx:start_idx + ns; 
    k = -Gz_amp./(ramp_samples(end) - ramp_samples(1));
    doff = Gz_amp - k*ramp_samples(1);
    Gz(ramp_samples) =k.*ramp_samples + doff;
    
    
    %Check if a polarity of the slice selection gradient is set. If not
    %default polarity is positive. (we need to invert previous calculated
    %Gz)
    if isfield(pulse, 'GsPolarity')
        switch pulse.GsPolarity
            case 'positive'
                Gz = (-1)*Gz;  
            case 'negative'
                Gz = Gz;  
            otherwise
                warning('pulse.GsPolarity can either be ''positive'' or ''negative''!');
        end
    else
        Gz = Gz*(-1);   
        warning('pulse.GsPolarity not set,  default is positive !');
    end

    
    %overlay Gmac gradient
    Gz = Gz + Gmac; 
    
    bplot = 0; 
    if bplot == 1
        figure; 
        subplot(2,1,1); 
        plot(d.tdis, B1*10);
        xlabel('time in ms'); 
        ylabel('B1 in �T'); 
        grid on
        title('RF Pulse'); 
        subplot(2,1,2); 
        plot(d.tdis, Gz); 
        if Gz_amp > 0 
            ylim([-1.4*Gz_amp, 1.4*Gz_amp]); 
        end
        grid on
        xlabel('time in ms'); 
        ylabel('Gz in mT/m'); 
        title('Slice selection gradient'); 
    end
    
    
    %% set parameters of Christoph Aigner's script 

    % space discretization
    d.a    = 0.10;                     % domain border in m
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

    d.u = B1(1,:); 
    d.v = B1(2,:); 
    
    % initial magnetization
    d.M0    = d.M0c*repmat([0;0;1],1,d.Nx); 

    M = cn_bloch(d,d.M0,d.u,d.v,d.w ); 

    Mcplx = squeeze(M(1:2,:,end)); 

    Mxy = sqrt(Mcplx(2,:).^2+Mcplx(1,:).^2); 
    
    try
      FWHM = fwhm(d.xdis*1E3, Mxy); 
    catch 
      FWHM = NaN; 
    end
    


    

end

