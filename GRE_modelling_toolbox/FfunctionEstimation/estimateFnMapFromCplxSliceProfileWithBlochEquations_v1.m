function [ Fn_mat] = estimateFnMapFromCplxSliceProfileWithBlochEquations_v1(pulse, Gsus, vx_size,  mask, te, bSliceScaling, B1_map)
%ESTIMATEFN Estimates the dephasing of spins by taking into account the
%slice profile, which is calculated by numeerical integration of the Bloch 
%equations, assuming a constant field gradient Gsus. 
%The simulation include the magnitude and the phase of the profile. It
%includes optional also B1 and the slice broading thinning due to the field
%gradient Gsus described by the scaling factor lamba = Gz_amp/(Gz_amp +
%Gsus). 
%The polarity of the slice select gradient is described withing the pulse
%struct with the field 'GsPolarity' in the pulse struct. 
%   @param   Gsus            Gradient map of B0 in z direction in mT/m [Ny, Nx, Nz]
%   @param   vx_size         Voxel size in z in mm 
%   @param   pulse           Struct with pulse parameters
%   @param   mask            Binary mask [Ny, Nx, Nz]
%   @param   te              Echo time in ms
%   @param   bSliceScaling   Boolean, if true lambda is calcualted
%   @param   B1-map          Normalized B1 map (optional)
%
% v1: Polarity is changed. Positive slice direction means that the it is
% positve in the gradient coordinate system, which is given by the
% cross-prodcut of the phase encdoing (e.g. R>>L) and the readout
% direction. Assuming that patient is head first the postive slice
% direciton is opposed to the device coordiante system where z points out
% of the scanner. 
%
%
%   Author: Martin Soellradl
%   Department of Neurology, Medical University of Graz, Graz, Austria
%   email:martin.soellradl@medunigraz.at
%   Website: http://www.neuroimaging.com
%   Januray 2020; Last revision: 02-Januray-2020

    % gyromagnetic ratio 42.57*2*pi in rad*1E6/T
    gamma = 267.51;                      
    

    %Check if B1-map is available
    if nargin < 7
       bB1corr = 0;  
    else
       bB1corr = 1;  
    end
    
    
    %Check if the influence of Gsus on the slice profile should be
    %included-
    %Apparently a negative polarity means that the gradient in z-direction
    %is positive and a positive polarity is negative in z-diretion. Weird. 
    if bSliceScaling == 1
        if isfield(pulse, 'GsPolarity')
            switch pulse.GsPolarity
                case 'positive'
                    polarity = 1;  
                case 'negative'
                    polarity = -1;  
                otherwise
                    error('pulse.GsPolarity can either be ''positive'' or ''negativ''!');
            end
        else
            %default polarity of slice select gradient is positive
            polarity = 1;  
        end
    else
        polarity = 1;   %needs to be defined for parfor
    end

    
    Nte = length(te); 
    [Ny, Nx, Nz] = size(mask); 
    
    %Calcualte slice selection gradient 
    Gz_amp = pulse.k_pulse/vx_size; 
    
    if bB1corr == 1
        %Simulate slice profile with Bloch equations for different B1 
        B1_steps = [0.5:0.05:1.5];
        %store nominal flip angle of the sequence
        alpha_nominal = pulse.alpha; 
        for i=1:length(B1_steps); 
            disp(['Bloch solver i is ' , num2str(i), ' out of ', num2str(length(B1_steps))]); 
            pulse.alpha = alpha_nominal.*B1_steps(i);  
            [Mcplx(i,:,:), ~, d] = simSliceProfileBlochEqn_v2(pulse, Gz_amp, 0);
        end
    else
        [Mcplx(1,:,:), ~, d] = simSliceProfileBlochEqn_v2(pulse, Gz_amp, 0);
    end
    
    %sampling points in mm 
    z = d.xdis*1E3; 
    
    %Crop profile 
    [~ ,idx_start] = min(abs(z + 4*vx_size)); 
    [~ ,idx_end] = min(abs(z - 4*vx_size)); 
    z_pulse = z(idx_start:idx_end); 
    Mcplx = Mcplx(:,:, idx_start:idx_end); 
    
    [theta, Mxy] = cart2pol(Mcplx(:,2,:), Mcplx(:,1,:));
    theta = squeeze(theta); 
    Mxy = squeeze(Mxy); 
    
    %get indices of non-zero voxels
    [idy,idx,idz] = ind2sub(size(mask),find(mask > 0));
    
    %number of non-zero voxles
    Nvx = length(idy); 
    
    k_lin = find(mask > 0);
    Gsus_vec = Gsus(k_lin); 
    
    if bB1corr == 1
        %Interpolate B1_steps
        X = d.xdis(idx_start:idx_end); 
        Y = B1_steps; 
        B1_intp = B1_steps(1):(B1_steps(end) -B1_steps(1))/100:B1_steps(end);
        [Xq, Yq] = meshgrid(X, B1_intp); 
        Mxy_intp = interp2(X,Y,Mxy,Xq,Yq);
        theta_intp = interp2(X,Y,theta,Xq,Yq);
        B1_vec = B1_map(k_lin); 
    else
        B1_vec = zeros(Nvx,1); B1_intp = 0; %init for parfor
        Mxy_intp = 0; theta_intp  = 0;
        idx_B1min = 0; 
        B1_vx = 0; 
        Mxy_mat = 0; 
        
    end

    %add te = 0 for normalization 
    te = [0, te]; 
    Nq = length(z_pulse); 
    Fn_vec = zeros(Nvx, Nte); 
    Fn_mat = zeros(Ny, Nx, Nz, Nte); 
    
    parfor vx=1:Nvx
        if ~mod(vx,5000)
            disp(['Voxel number is: ', num2str(vx), ' out of: ', num2str(Nvx)]); 
        end

        Gsus_vx = Gsus_vec(vx); %mT/m
        
        %Select the best slice profile that fits B1 value
        if bB1corr == 1
            B1_vx = B1_vec(vx); 
            [m, idx_B1min] = min(abs(B1_intp - B1_vx)); 
            Mxy_mat = repmat(Mxy_intp(idx_B1min, :), [Nte + 1, 1]); 
            theta_mat = repmat(theta_intp(idx_B1min, :), [Nte + 1, 1]); 
        else
            Mxy_mat = repmat(Mxy', [Nte + 1, 1]);
            theta_mat = repmat(theta', [Nte + 1, 1]); 
        end
        
        %adapt the slice width based on the field gradient 
        z = z_pulse; 
        if bSliceScaling == 1
           lambda = polarity.*Gz_amp./(polarity.*Gz_amp + Gsus_vx); 
           z = z_pulse.*lambda; 
        end
 
        wz = gamma*Gsus_vx.*z; % [2*pi MHz/T * mT/m *mm] = rad/s
        wz_mat = repmat(wz, [Nte + 1, 1]); %Nte + 1 because of normalization
        te_mat = repmat(te', [1, Nq]); 

        %dephasing weighted with the profile and integrate
        Mxy_te = Mxy_mat.*exp(1i*wz_mat.*te_mat*1E-3 + 1i.*theta_mat);
        S = squeeze(sum(Mxy_te, 2)); 

        %normalize signal
        S0 = abs(S(1)); 
        Snorm = S./S0; 
        Fn_vec(vx, :) = abs(Snorm(2:end));   

    end
    
    for vx=1:Nvx
        Fn_mat(idy(vx), idx(vx),  idz(vx), :) = squeeze(Fn_vec(vx, :)) ;   
    end


end

