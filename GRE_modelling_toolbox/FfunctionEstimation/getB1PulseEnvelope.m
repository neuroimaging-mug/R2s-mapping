function [ B1 ] = getB1PulseEnvelope( pulse, t, Tramp)
%GETB1PULSEENVELOPE Returns Pulse Envelope for given pulse
%   @param      pulse   Struct with pulse parameters
%   @param      t       Sampling points in ms 
%   @param      Tramp   Ramp-Up/Down time of gradient (optional)

%   Author: Martin Soellradl
%   Department of Neurology, Medical University of Graz, Graz, Austria
%   email:martin.soellradl@medunigraz.at
%   Website: http://www.neuroimaging.com
%   Januray 2020; Last revision: 02-Januray-2020

    if nargin < 3
        Tramp = 0; 
    end
    
    gamma = 267.51;  % 42.57*2*pi in MHz/T
    
    ts = t(2) - t(1); %sampling interval
    %create a window where RF-pulse is switched on
    pulse_wdw = zeros(1, length(t)); 
    pulse_wdw(round(Tramp/ts): round(Tramp/ts)+ round(pulse.Tpulse/ts)) = 1; 
    
    %time offset of pulse and window
    tshift = (t  - pulse.Tpulse/2 - Tramp); 
    B1 = zeros(2, length(t)); % B1(1,:) real part, B1(2,:) imag
    switch pulse.type
       case 'sinc-hanning'
            
            wdw = (1 + cos(2*pi*(tshift /pulse.Tpulse)))*0.5.*pulse_wdw; %hanning window times pulse window
            %Attentention with sinc definition of maltab sinc(t) = sin(pi*t)/(pi*t) for t not 0
            sinc_fun = sinc(pulse.BWT.*tshift/pulse.Tpulse);
            %sinc-hanning function without flip angle scaling
            B1(1,:) = sinc_fun.*wdw; 
            %estimate correct scaling to get desired flip angle
            Amp =  pulse.alpha/rad2deg(gamma.*trapz(t*1E-3,   squeeze(B1(1,:))*10)); %weird unit thing [B1] = 0.1 �T
            B1(1,:) = Amp*B1(1,:);
        case 'hanning-wdw'
            wdw = (1 + cos(2*pi*(tshift /pulse.Tpulse)))*0.5.*pulse_wdw; %hanning window times pulse window
            %Attentention with sinc definition of maltab sinc(t) = sin(pi*t)/(pi*t) for t not 0
            %sinc-hanning function without flip angle scaling
            B1(1,:) = wdw; 
            %estimate correct scaling to get desired flip angle
            Amp =  pulse.alpha/rad2deg(gamma.*trapz(t*1E-3,   squeeze(B1(1,:))*10)); %weird unit thing [B1] = 0.1 �T
            B1(1,:) = Amp*B1(1,:);
       case 'gauss'
            %gaussian without flip angle scaling
            B1(1,:) = exp(-tshift.^2./(2*pulse.sigma .^2)); 
            %estimate correct scaling to get desired flip angle
            Amp =  pulse.alpha/rad2deg(gamma.*trapz(t*1E-3,   B1(1,:)*10)); %weird unit thing [B1] = 0.1 �T
            B1(1,:) = Amp*B1(1,:).*pulse_wdw;
            
        case 'exponential'
            B1(1,:) = exp(- 8*abs(tshift)/pulse.Tpulse); 
            Amp =  pulse.alpha/rad2deg(gamma.*trapz(t*1E-3,   B1(1,:)*10)); %weird unit thing [B1] = 0.1 �T
            B1(1,:) = Amp*B1(1,:).*pulse_wdw; 
        case 'NL_filter'
            wdw = (1 + cos(2*pi*(tshift /pulse.Tpulse)))*0.5; 
            
            B1_mag =  sin(2*pi.*tshift./pulse.Tpulse).*pulse_wdw.*wdw;
            [B1(1,:), B1(2,:)] =  pol2cart(pi/2, B1_mag); 
            B1(2,:) = B1(2,:)*(-1); 
            %estimate correct scaling to get desired flip angle
            Amp =  pulse.alpha/rad2deg(gamma.*trapz(t*1E-3, abs(B1(2,:))*10)); %weird unit thing [B1] = 0.1 �T
            B1(2,:) = Amp*B1(2,:).*pulse_wdw;
            
        otherwise
          disp('Unkown pulse type'); 
    end
end

