%Script plots the slice profile for three different sinc-pulses with the
%phase. 
% If you use this code please cite: 
%
%   Soellradl M, Lesch A, Strasser J, et al. 
%   Assessment and correction of macroscopic field variations in 2D spoiled 
%   gradient-echo sequences. Magn Reson Med. 2019;00:114
%   https ://doi.org/10.1002/mrm.28139
%
% Author: Martin Soellradl
% Department of Neurology, Medical University of Graz, Graz, Austria
% email:martin.soellradl@medunigraz.at
% Website: http://www.neuroimaging.com
% Januray 2020; Last revision: 02-Januray-2020



addpath('external_toolboxes'); 
addpath('GRE_modelling_toolbox'); 

alpha = [30,60, 90]; %deg 
vx_size = 4;         %mm


%% Load RF-fast (BWT=2, Tpulse =1ms)
    
pulse.type = 'sinc-hanning'; 
pulse.Tpulse = 1;               %ms 
pulse.BWT = 2;                %-
pulse.k_pulse = 56.6132;         %�T
pulse.name  ='sinc-hanning BWT = 2 Tpulse = 1ms'; 

pulse_cell{1} = pulse; 
 
%% RF normal (sinc-hanning BWT=2.7, Tpulse = 2ms)

pulse.type = 'sinc-hanning'; 
pulse.Tpulse = 2;               %ms 
pulse.BWT = 2.7;                %-
pulse.k_pulse = 33.1510;        %�T/m
pulse.name  ='sinc-hanning BWT = 2.7 Tpulse = 2ms'; 
pulse_cell{2} = pulse; 

%% RF long (sinc-hanning BWT=8, Tpulse = 4ms) 

pulse.type = 'sinc-hanning'; 
pulse.Tpulse = 4;               %ms 
pulse.BWT = 8;                %-
pulse.k_pulse = 46.9503;       %�T/m
pulse.name  ='sinc-hanning BWT = 8 Tpulse = 4ms'; 

pulse_cell{3} = pulse; 




%% Do the main simulation 

polarity = {'positive', 'negative'}; 

for k=1:2; 
    for i=1 %:length(pulse_cell); 
       disp(['Pulse nr: ' , num2str(i), ' out of ' num2str(length(pulse_cell))]); 
       for j= 1:length(alpha); 
           pulse_tmp = pulse_cell{i}; 
           pulse_tmp.alpha = alpha(j); 
           pulse_tmp.GsPolarity = polarity{k}; 
           Gz_amp = pulse_tmp.k_pulse/vx_size; 

           [Mcplx(k,i,j,:,:), ~, d] = simSliceProfileBlochEqn_v2(pulse_tmp, Gz_amp, 0);
       end
    end
end


%%

%sampling points in mm 
z = d.xdis*1E3; 

%Crop profile 
[~ ,idx_start] = min(abs(z + vx_size)); 
[~ ,idx_end] = min(abs(z - vx_size)); 
z_pulse = z(idx_start:idx_end); 
Mcplx_crp = Mcplx(:,:,:,:,idx_start:idx_end); 

[theta, Mxy] = cart2pol(Mcplx_crp(:,:,:,2,:), Mcplx_crp(:,:,:,1,:));
theta = squeeze(theta); 
Mxy = squeeze(Mxy); 


%% store results

path = 'results/PulseComparisionWithPhase/',
if exist(path, 'dir') == 0
            mkdir(path);   
end


%store results 
save([path, 'Mcplx_crp'], 'Mcplx_crp'); 
save([path, 'theta'], 'theta'); 
save([path, 'Mxy'], 'Mxy'); 
save([path, 'z_pulse'], 'z_pulse'); 
save([path, 'alpha'], 'alpha'); 
save([path, 'pulse_cell'], 'pulse_cell'); 

%% Plot results


Npulse = length(pulse_cell); 
polrity_correct = {'positive', 'negative'}; 

for k=1:2
    fig = figure; 
    screensize = get( 0, 'Screensize' );
%     set(fig,'Position',screensize);
%     set(fig, 'PaperPositionMode', 'auto')
    for i=1:Npulse
      pulse_tmp = pulse_cell{i}; 
      subplot(3,2,2*i -1); 
        plot(z_pulse, squeeze(Mxy(k,i,:,:)), 'LineWidth', 2); 
        xlabel('z pos in mm'); ylabel('|M_{xy}(z)|'); 
        title(['Magnitude ' , pulse_tmp.name, ' Gs ', polrity_correct{k}]);
      subplot(3,2,2*i); 
        plot(z_pulse, squeeze(theta(k,i,:,:)), 'LineWidth', 2); 
        xlabel('z pos in mm'); ylabel('\phi_{xy}(z)'); 
        title(['Phase ' , pulse_tmp.name, ' Gs ', polrity_correct{k}]);
        xlim([-4,4]);  ylim([-0.4, 0.4]); 
    end

    fig_name = [path, 'PulseComparision_Gs_', polrity_correct{k}];  
%     export_fig(fig_name, '-png'); 

    print(fig_name,'-dpng',  '-r300'); % print(fig_name,'-depsc');
end


