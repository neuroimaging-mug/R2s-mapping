function [R2s_corr, A, resnorm_map, residuals_map] = CorrectedR2sMapEstimationFnNonLinFit(mag, Fn, te, mask)
%R2SMAPESTIMATION Estimates correted R2smap from mgre data.
%   Script performs an corrected R2s estimation with lsqnonlin.
%   It estimates R2s and A with the following the following model from 
%   Smodel= A*exp(-t*R2s)*Fn
%   @param      mag     Magnitude image [Nx, Ny, Nz, t]
%   @param      Fn      F-function describes the
%                       signal behavior due to macroscopic field 
%                       inhomogeneities 
%   @param      te      echo times in ms 
%   @param      mask    Binary mask with [Nx, Ny, Nz]
%
%   Author: Martin Soellradl
%   Department of Neurology, Medical University of Graz, Graz, Austria
%   email:martin.soellradl@medunigraz.at
%   Website: http://www.neuroimaging.com
%   November 2019; Last revision: 01-November-2019

   
    [idy,idx,idz] = ind2sub(size(mask),find(mask > 0));
    
    N = size(mask); 
    
    if length(N) == 2
        N(3) = 1;  
        mag = permute(mag, [1,2,4, 3]); 
        Fn = permute(Fn, [1,2,4, 3]); 
    end
    
    %init result matrices 
    fit_results = zeros([N, 2]); 
    resnorm_map  = zeros(N); 
    residuals_map = zeros([N, length(te)]); 
    fit_solution = zeros(length(idz),2); 
    resnorm = zeros(length(idz),1); 
    residual = zeros(length(idz),length(te)); 
    fit_prop.name = 'Exponential Fit';  

    Nvx = length(idz);
    parfor vx=1:Nvx; 
         if ~mod(vx,5000)
              disp(['Voxel number is: ', num2str(vx), ' out of: ', num2str(Nvx)]); 
         end
         Smeas = squeeze(mag(idy(vx), idx(vx),idz(vx), :))'; 
         [fit_solution(vx,:), resnorm(vx),residual(vx,:)] ...
              =  ExponentialR2starFitWithFnFunction(te, Smeas, squeeze(Fn(idy(vx), idx(vx),idz(vx),:)), fit_prop,0);

    end

     %write fit results into image 
    for vx = 1:length(idz); 
        fit_results(idy(vx),idx(vx), idz(vx),:) = fit_solution(vx,:); 
        residuals_map(idy(vx),idx(vx), idz(vx),:) = residual(vx,:);
        resnorm_map(idy(vx),idx(vx), idz(vx)) = resnorm(vx);
    end
    
    A = squeeze(fit_results(:,:,:,1)); 
    R2s_corr = squeeze(fit_results(:,:,:,2));
    
    


end

