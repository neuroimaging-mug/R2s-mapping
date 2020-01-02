%Example 1 demonstrates the influence of the slice selection polarity Gslice 
%in presence of a macrospoic field gradient Gz. 

pulse.type = 'sinc-hanning'; 
pulse.Tpulse = 2;               %ms 
pulse.BWT = 2.7;                %-
pulse.k_pulse = 33.1510;        %�T/m
pulse.name  ='sinc-hanning BWT = 2.7 Tpulse = 2ms'; 
z0 = 4; %mm
alpha = [30, 90]; %deg 
pulse_cell{1} = pulse; 

polarity = {'positive', 'negative'}; 
for k=1:2; 
    for i=1 :length(pulse_cell); 
       disp(['Pulse nr: ' , num2str(i), ' out of ' num2str(length(pulse_cell))]); 
       for j= 1:length(alpha); 
           pulse_tmp = pulse_cell{i}; 
           pulse_tmp.alpha = alpha(j); 
           pulse_tmp.GsPolarity = polarity{k}; 
           Gz_amp = pulse_tmp.k_pulse/z0; 

           [Mcplx(k,i,j,:,:), ~, d{k,i,j}] = simSliceProfileBlochEqn_v2(pulse_tmp, Gz_amp, 0);
       end
    end
end


%sampling points in mm 
z = d{1,1,1}.xdis*1E3; 

vx_size = z0; 
%Crop profile 
[~ ,idx_start] = min(abs(z + vx_size)); 
[~ ,idx_end] = min(abs(z - vx_size)); 
z_pulse = z(idx_start:idx_end); 
Mcplx_crp = Mcplx(:,:,:,:,idx_start:idx_end); 

[theta, Mxy] = cart2pol(Mcplx_crp(:,:,:,2,:), Mcplx_crp(:,:,:,1,:));
theta = squeeze(theta); 
Mxy = squeeze(Mxy); 


path = 'Simulations/mat_files/',
if exist(path, 'dir') == 0
            mkdir(path);   
end


%store results 
%store results 
save([path, 'd'], 'd'); 
save([path, 'Mcplx'], 'Mcplx'); 
save([path, 'Mcplx_crp'], 'Mcplx_crp'); 
save([path, 'theta'], 'theta'); 
save([path, 'Mxy'], 'Mxy'); 
save([path, 'z_pulse'], 'z_pulse'); 
save([path, 'alpha'], 'alpha'); 
save([path, 'pulse_cell'], 'pulse_cell'); 

%% Sim dephasing 
Gsus_vx = 0.1 %mT/m
gamma = 267.51;                      % gyromagnetic ratio 42.57*2*pi in MHz/T
te = [0:0.1:100]; %ms
Nte = length(te); 
Nq = length(z_pulse); 

for k=1:2; 
        for j= 1:length(alpha); 
            Mxy_sim = squeeze(Mxy(k,j,:,:)); 
            Mcplx_mat = repmat(Mxy_sim', [Nte, 1]);
            wz = gamma*Gsus_vx.*z_pulse; % [2*pi MHz/T * mT/m *mm] = rad/s
            wz_mat = repmat(wz, [Nte, 1]); %Nte + 1 because of normalization
            te_mat = repmat(te', [1, Nq]); 

            theta_mat = squeeze(repmat(theta(k,j,:,:), [Nte, 1])); 
            %dephasing weighted with the profile and integrate
            Mxy_te = Mcplx_mat.*exp(1i*wz_mat.*te_mat*1E-3 + 1i.*theta_mat);
            S = squeeze(sum(Mxy_te, 2)); 

            %normalize signal
            S0 = abs(S(1)); 
            Snorm(k,j,:) = S./S0; 
        end
end


%% 

path = 'Simulations/Figures/PulseComparisionWithPhase/',
if exist(path, 'dir') == 0
            mkdir(path);   
end

F_size = 20; 
fig = figure;
screensize = get( groot, 'Screensize' );
set(fig,'Position',screensize);
set(fig, 'PaperPositionMode', 'auto')
subplot(2,3,1); 
    plot(z_pulse, squeeze(Mxy(1,1,:,:)), 'LineWidth', 3); hold on; 
    plot(z_pulse, squeeze(Mxy(2,1,:,:)), '--',  'LineWidth', 4); 
    xlabel('z pos (mm)'); ylabel('|M_{xy}| (a.u.)'); 
    title('Magnitude, \alpha = 30�');xlim([-4,4]);
    legend('G_{slice} neg.', 'G_{slice} pos.');
    grid on; 
subplot(2,3,2); 
    plot(z_pulse,squeeze(theta(1, 1,:)), 'LineWidth', 4); hold on; 
    plot(z_pulse,squeeze(theta(2, 1,:)), '--', 'LineWidth', 3); 
    xlabel('z pos (mm)'); ylabel('\phi_{xy} (rad)'); 
    title('Phase, \alpha = 30�');
    xlim([-4,4]); ylim([-pi/8, pi/8]); 
    legend('G_{slice} neg.', 'G_{slice} pos.');
    grid on; 
subplot(2,3,4); 
    plot(z_pulse, squeeze(Mxy(1,2,:,:)), 'LineWidth', 3); hold on; 
    plot(z_pulse, squeeze(Mxy(2,2,:,:)), '--',  'LineWidth', 3); 
    xlabel('z pos (mm)'); ylabel('|M_{xy}| (a.u.)'); 
    title('Magnitude, \alpha = 90�');xlim([-4,4]);
    legend('G_{slice} neg.', 'G_{slice} pos.');
    grid on; 
subplot(2,3,5); 
    plot(z_pulse,squeeze(theta(1,2,:)), 'LineWidth', 3); hold on; 
    plot(z_pulse,squeeze(theta(2,2,:)), '--', 'LineWidth', 3); 
    
    xlabel('z pos (mm)'); ylabel('\phi_{xy} (rad)'); 
    title('Phase, \alpha = 90�');
    xlim([-4,4]);  ylim([-pi/8, pi/8]); 
    legend('G_{slice} neg.', 'G_{slice} pos.');
    grid on;  
subplot(2,3,3); 
    plot(te, abs(squeeze(Snorm(1,1,:,:))), 'LineWidth', 3); hold on; 
    plot(te, abs(squeeze(Snorm(2,1,:,:))),'--', 'LineWidth', 3); 
    xlabel('TE (ms)'); ylabel('F (a.u.)'); 
    title(['Gz = ', num2str(Gsus_vx*1E3), '�T/m', ' \alpha = 30�']);
    legend('G_{slice} neg.', 'G_{slice} pos.');
    grid on; 
    ylim([0,1.1]); 
subplot(2,3,6); 
    plot(te, abs(squeeze(Snorm(1,2,:,:))), 'LineWidth', 3); hold on; 
    plot(te, abs(squeeze(Snorm(2,2,:,:))),'--', 'LineWidth', 3); 
    xlabel('TE (ms)'); ylabel('F (a.u.)'); 
    title(['Gz = ', num2str(Gsus_vx*1E3), '�T/m', ' \alpha = 90�']);
    legend('G_{slice} neg.', 'G_{slice} pos.');
    grid on; 
        ylim([0,1.1]); 
    
set(findall(gcf,'-property','FontSize'),'FontSize',F_size)
 
fig_name = [path, 'sinc_normal_30_90deg']; 
fig_name = [path, 'sinc_normal_30_90deg_print']; 
print(fig_name,'-dpng',  '-r300');  print(fig_name,'-depsc'); 
    print(fig_name,'-dsvg');

