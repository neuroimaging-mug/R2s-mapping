function [ Fn_wSlcCorr, Fn_preibisch, B1_env ] = estimateFnMapWtihRFShape( pulse, Gsus, te, z0 )
%ESTIMATEFNMAPWTIHRFSHAPE Estimates the signal attenuation of the gradient
%echo signal for a macroscopic field gradient Gsus depending on the
%RF-pulse shape and the echo times as described in Preibisch et al MRM
%2008.
%   @param      pulse   Struct with pulse parameters
%   @param      Gsus    Macroscopic field gradient in z direction [mT/m] 
%   @param      te      Echo times in ms 
%   @param      z0      Slice thickness in mm
%
%
%   Author: Martin Soellradl
%   Department of Neurology, Medical University of Graz, Graz, Austria
%   email:martin.soellradl@medunigraz.at
%   Website: http://www.neuroimaging.com
%   Januray 2019; Last revision: 01-Januray-2019


    %get pulse functin and the specific pulse constant for the calcuation
    %of the slice selection gradient (depends on pulse properties!)
    switch pulse.type
       case 'sinc-hanning'
            B1 = @(t) sincHanningRecWdw(t, pulse.BWT, pulse.Tpulse); 
            
            if pulse.BWT == 2.7 && pulse.Tpulse == 2
               k_pulse = 33.1510;
            elseif  pulse.BWT == 2 && pulse.Tpulse == 1
               k_pulse = 56.6132;
            elseif  pulse.BWT == 8 && pulse.Tpulse == 4
               k_pulse = 46.9503;
            end
            
       case 'gauss'
            B1 = @(t) exp(-t.^2./(2*pulse.sigma .^2)); 
            
            if pulse.Tpulse == 2 && pulse.sigma == 0.280
               k_pulse = 31.6927; 
            elseif pulse.Tpulse == 1  pulse.sigma = 0.162
               k_pulse = 54.9298; 
            end
            
       case 'exponential'
            B1 = @(t) exp(-8*abs(t)./pulse.Tpulse); 
            
            k_pulse = 30.6553; 
       otherwise
          disp('Unkown pulse type'); 
    end
    
    Gslice = k_pulse/z0; %mT/m
    
    [Ny, Nx, Nz] = size(Gsus); 
    Nte = length(te); 
    
    te_mat = repmat(permute(te(:), [4,2,3,1]), [Ny, Nx, Nz, 1]); 
    Gsus_mat = repmat(Gsus, [1,1,1,Nte]); 
    Gslice_mat = repmat(Gslice, [Ny,Nx, Nz, Nte]); 
    
    tsubs_wSlcCorr = (Gsus_mat./(Gslice_mat - Gsus_mat)).*te_mat;
    tsubs_preibisch = (Gsus_mat./Gslice_mat).*te_mat;
    
%     figure; 
%     subplot(1,2,1); plot(squeeze(tsubs_wSlcCorr), squeeze(B1(tsubs_wSlcCorr))); 
%     subplot(1,2,2); plot(abs(fftshift(fft(squeeze(B1(tsubs_wSlcCorr)))))); 



    Fn_wSlcCorr = B1(tsubs_wSlcCorr); 
    Fn_preibisch = B1(tsubs_preibisch); 
    B1_env = B1(te);

end


function B1 = sincHanningRecWdw(t, BWT, Tpulse)

    rec_wdw = ones(size(t));  
    rec_wdw(abs(t) > Tpulse/2) = 0;
    
    B1 = (1 + cos(2*pi*(t./Tpulse)))*0.5...
                .*sinc(BWT.*t/Tpulse).*rec_wdw; 
end
