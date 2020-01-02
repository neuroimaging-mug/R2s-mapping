function [ dw0, phi0, resnorm ] = CplxFitOfB0map(phase,mask, te_selected, bplot)
%LINEARFITOFB0MAP Conplex fit of the phase signal to estimate dw0
%   @param      phase            Phaese image [Nx, Ny, Nz, t]
%   @param      te_selected      echo times [ms] 
%   @param      mask             [Nx, Ny, Nz]
%   @param      bplot            Display fit result for voxel
%   @param      dw0              Frequncy map [rad/s]
%   @param      phi0             Phase offset [rad]
%
%
%   Author: Martin Soellradl
%   Department of Neurology, Medical University of Graz, Graz, Austria
%   email:martin.soellradl@medunigraz.at
%   Website: http://www.neuroimaging.com
%   Januray 2019; Last revision: 01-Januray-2019

    te_selected = te_selected(:); 
    te = te_selected; 
    t_fit = 0:te(end)/100:te(end); 
    %field assumption 
    x0 = [0,0]; opts = optimset('Display','off');
    fun = @(x,te_val) exp(1j*(x(1) + x(2).*te_val)); 
    lb = [-pi -1]; 
    ub = [pi, 1]; 

    sel_slc = [1:size(phase,3)]; 

    %initialization of coeffcient map 
    par_est = zeros(size(phase,1), size(phase,2), length(sel_slc),2); 
    resnorm = zeros(size(phase,1), size(phase,2), length(sel_slc)); 
    for i= 1:length(sel_slc); 
        disp(['Slice number is: ', num2str(i), ' out of: ', num2str(length(sel_slc))]); 
        img_phase_slc = squeeze(phase(:,:,sel_slc(i),:));
        [idy, idx] = find(mask(:,:,i) > 0);
        par_tmp = zeros(length(idx),2); 
        resnorm_tmp = zeros(length(idx),1); 
        parfor vx = 1:length(idx)
            phase_time = exp(1j*squeeze(img_phase_slc(idy(vx),idx(vx),1:length(te)))); 
            %estimate parameters of line
            [x_est, resnorm_tmp(vx,:)] = lsqcurvefit(fun,x0,te, phase_time ,lb,ub,opts);
            par_tmp(vx,:) = x_est;  
            if bplot
                t_fit = 0:te(end)/50:te(end); 
                figure; 
                scatter(te,angle(phase_time));hold on; 
                title('Linear Fit Approach'); 
                plot(t_fit, angle(fun(x_est, t_fit))); 
            end
        end
        for vx=1:length(idx)
            par_est(idy(vx), idx(vx),i,:) = par_tmp(vx,:); 
            resnorm(idy(vx), idx(vx),i) =  resnorm_tmp(vx,:);
        end
    end
    phi0 = par_est(:,:,:,1); 
    dw0 = par_est(:,:,:,2)*1E3; %rad/s
end

