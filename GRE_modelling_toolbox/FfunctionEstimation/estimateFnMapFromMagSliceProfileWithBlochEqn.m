function [ Fn_mat] = estimateFnMapFromMagSliceProfileWithBlochEqn(pulse, Gsus, vx_size,  mask, te, bSliceScaling)
%ESTIMATEFN Estimates the dephasing of spins by taking into account the
%slice profile, which is calculated by numeerical integration of the Bloch 
%equations, assuming a constant field gradient Gsus. 
%   @param   Gsus            Gradient map of B0 in z direction in mT/m [Ny, Nx, Nz]
%   @param   vx_size         Voxel size in z 
%   @param   pulse           Struct with pulse parameters
%   @param   mask            Binary mask [Ny, Nx, Nz]
%   @param   te              Echo time in ms


%   Author: Martin Soellradl
%   Department of Neurology, Medical University of Graz, Graz, Austria
%   email:martin.soellradl@medunigraz.at
%   Website: http://www.neuroimaging.com
%   November 2019; Last revision: 01-November-2019
  
    if nargin < 6
        bSliceScaling = 0; 
    end
    
    % gyromagnetic ratio 42.57*2*pi in MHz/T
    gamma = 267.51;                      
    
    Nte = length(te); 
    [Ny, Nx, Nz] = size(mask); 
    
    %Calcualte slice selection gradient 
    Gz_amp = pulse.k_pulse/vx_size; 
    
    %Simulate slice profile with Bloch equations
    if bSliceScaling == 1
        [Mcplx, ~, d] = simSliceProfileBlochEqn_v2(pulse, Gz_amp, Gsus);  
    else
        [Mcplx, ~, d] = simSliceProfileBlochEqn_v2(pulse, Gz_amp, 0);  
    end
    
    %sampling points in mm 
    z = d.xdis*1E3; 

    %Crop profile 
    [~ ,idx_start] = min(abs(z + 8*vx_size)); 
    [~ ,idx_end] = min(abs(z - 8*vx_size)); 
    z_pulse = z(idx_start:idx_end); 
    Mcplx = Mcplx(:,idx_start:idx_end); 
    
    %use Mx (real part) of the slice profile to calculate dephasing
%     Mxy = Mcplx(2,:); 
    Mxy = sqrt(Mcplx(2,:).^2+Mcplx(1,:).^2); 
    
    [idy,idx,idz] = ind2sub(size(mask),find(mask > 0));
    
    k_lin = find(mask > 0);
    Gsus_vec = Gsus(k_lin); 
    
    Nvx = length(idy); 
    
    %add te = 0 for normalization 
    te = [0, te]; 
    
    Mxy_mat = repmat(Mxy, [Nte + 1, 1]); 
    Nq = length(z_pulse); 
    Fn_vec = zeros(Nvx, Nte); 
    Fn_mat = zeros(Ny, Nx, Nz, Nte); 
    

    
    parfor vx=1:Nvx
        if ~mod(vx,5000)
            disp(['Voxel number is: ', num2str(vx), ' out of: ', num2str(Nvx)]); 
        end

        Gsus_vx = Gsus_vec(vx); %mT/m

        wz = gamma*Gsus_vx.*z_pulse % [2*pi MHz/T * mT/m *mm] = rad/s

        wz_mat = repmat(wz, [Nte + 1, 1]); %Nte + 1 because of normalization
        te_mat = repmat(te', [1, Nq]); 

        %dephasing weighted with the profile 
        Mxy_te = Mxy_mat.*exp(-1i*wz_mat.*te_mat*1E-3);

        %integrate signal along z
        S = squeeze(trapz(Mxy_te, 2)); 

        %normalize signal
        S0 = abs(S(1)); 
        Snorm = S./S0; 

        Fn_vec(vx, :) = abs(Snorm(2:end));   

    end
    
    for vx=1:Nvx
        if ~mod(vx,5000)
            disp(['Writing into matrix.. voxel number is: ', num2str(vx), ' out of: ', num2str(Nvx)]); 
        end
        Fn_mat(idy(vx), idx(vx),  idz(vx), :) = squeeze(Fn_vec(vx, :)) ;   
    end


end

