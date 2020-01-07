function [MWF_map, M0_map, mask_signal] = NNLS_MultiExponentialMapEstimation_MERA_withCorrection(mag, te, mask, Fn, e_inv_th, beta_reg)
%R2SMAPESTIMATION Estimates MWF estimation with 
%a non-negatvie least squares (NNLS) multi-compartment

%   Script that performs an R2s estimation from multi gradient echo data. 
%   @param      mag       Magnitude image [Nx, Ny, Nz, t]
%   @param      te        echo times in ms 
%   @param      mask    
%   @param      Fn        Signal dephasing due to macroscopic field gradient Gz[Nx, Ny, Nz, t]
%   @param      e_inv_th  Error threshold for dividing Smeas./Fn

% Author: Martin Soellradl
% Department of Neurology, Medical University of Graz, Graz, Austria
% email:martin.soellradl@medunigraz.at
% Website: http://www.neuroimaging.com
% Januray 2020; Last revision: 02-Januray-2020


    th_MWF = 0.025; %s
    te = te*1E-3;  %must be in seconds! 
    
    if nargin < 6
        beta_reg = 0.001; 
    end
    
    fitting.regtyp = 'me';
    % % optional inputs for 'mg' fitting
    % fitting.widthgauss = 5;
    % fitting.numbergauss = 2;
    % fitting.T0gauss = [0.02 0.07]';
    % fitting.numberT0 = 5;

    % for conventional regularization, the regulization weighting must be
    % adjusted. This can be done manually (fitting.regadj = 'manual') or using
    % one of two automatic methods ('gcv' and 'erinc') -- see literature
    % references in code
    fitting.regadj = 'manual';
%     fitting.regweight = 0.001;
    fitting.regweight = beta_reg; %0.000001;

    % graph the results or not. This input is irrelevant if
    % analysis.interactive = 'y'

    analysis.graph = 'n';

    % define the range of time constants to fit. Note that for 'mg' fitting,
    % this is the full range of the spectrum, but the lowest and highest mean
    % component time constants (echoed to display) cover a narrower domain
    fitting.rangeT = [2.5e-3 0.25];

    % define the number of time constants in the spectrum
    fitting.numberT = 200;

    % set the non-negative least square (NNLS) code. In most cases, the supplied
    % nnlsmex fuction is the fastest NNLS solver. You may need to compile this
    % for your system. Or you can specify the MATLAB fuction,
    % fitting.nnlscode = 'lsqnonneg'. 

    % You can automatically or manually extract a finite number of component
    % amplitudes and time constants. Test this out with the interactive GUI.
    analysis.extract = 'auto';
    analysis.numberextract = 2;


    
    
    
    
    
    
    mask_noise = zeros(size(mask));
    mask_noise(1:20,1:20, :) = 1; 
    for i=1:size(mag,4);
      mag_tmp =mag(:,:,:,i); 
      std_noise(i) = std(mag_tmp(mask_noise == 1)); 
    end


    [Nx, Ny, Nz, NE] = size(mag);

    MWF_map = zeros(Nx, Ny, Nz); 
    M0_map = MWF_map; 
    mask_signal = zeros(Nx, Ny, Nz, NE); 
    SNR_min = 50; 
%     e_inv_th = e_inv_th; %0.05; 
    for i= 1:Nz 
       disp(['slice i: ', num2str(i), ' out of', num2str(Nz)]); 

       mag_slc =squeeze(mag(:,:,i,:)); 
       Fn_slc = squeeze(Fn(:,:,i,:)); 
       mask_slc = squeeze(mask(:,:,i)); 

       %Estimate SNR in each voxel 
       SNR_img = mag_slc(:,:,1)./std_noise(1); 

       %mask for inversion: 1) enough SNR 2) within the mask 
       SNR_mask = SNR_img; SNR_mask(SNR_mask < SNR_min) = 0; 
       SNR_mask(SNR_mask >= SNR_min) = 1;  
       SNR_mask = SNR_mask.*mask_slc; 

       %Do an error propagation
       d_Smeas = std_noise(1)./Fn_slc;  
       mag_corr =  mag_slc./Fn_slc; 
       mag_corr = mag_corr.*repmat(SNR_mask, [1,1, NE]); 
       e_inv = abs(d_Smeas)./mag_corr;

  

       %Now estimate for each voxel the echo time when error is biggher than
       %the error threshold. 
       e_inv_cum = cumsum(e_inv,3);

       %generate map with number of reliable echoes  
       mask_signal_slc = thresholdImage(e_inv_cum, [0, e_inv_th]); %2
       mask_signal(:,:,i,:) = mask_signal_slc;
       mask_echo = sum(mask_signal_slc, 3); 

       
       %Now do an NNLS fit for all realiable numbers of echoes
       nTE_low = 6; 

       mag_c_re = reshape(mag_corr, [Nx*Ny, NE]); 


       %seperate data for differnt echoes and prepare for parfor 
       for n=nTE_low:NE

         %numbers of echoes used for this iteration
         data.t = te(1:n)';

         %get indices that contain n echos
         k = find(mask_echo == n); 
         mag_tmp= mag_c_re(k, 1:n);
         data.D = mag_tmp'; 

         data_c{n} = data; 
         k_c{n} = k; 
       end


       %Now parfor (not really elegant...) 
       M0_slc = zeros(Nx, Ny); %sliced variables
       MWF_slc = M0_slc; 
       parfor n=nTE_low:NE %activate parfor

         data_tmp =  data_c{n}; %current data
         if ~isempty(data_tmp.D); 
             tic
                [out1D]=MERA(data_tmp,fitting,analysis);
             toc

             %The fitted signal amplitude at t = 0
             M0_vec = sum(out1D.S);

             %Estimate MWF 
             [~, idx_th] = min(abs( out1D.T - th_MWF)); 
             MWF_vec = sum(out1D.S(1:idx_th,:))./sum(out1D.S);

             MWF = zeros(Nx, Ny);   
             M0 = MWF; 
             k = k_c{n}; 
             for vx=1:length(k); 
                MWF(k(vx)) = MWF_vec(vx); 
                M0(k(vx)) =  M0_vec(vx); 
             end

             MWF_slc = MWF + MWF_slc;
             M0_slc = M0_slc + M0; 
         end 

       end

       M0_map(:,:,i) = M0_slc; 
       MWF_map(:,:,i) = MWF_slc; 


    end

   
    
    
    
   
end

