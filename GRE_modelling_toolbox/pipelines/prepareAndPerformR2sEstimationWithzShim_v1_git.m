function [ data ] = prepareAndPerformR2sEstimationWithzShim_v1_git( par )
%PREPAREANDPERFORMMWFESTIMATION Script loads raw data of the multi-echo
%gradient-echo data with z-shim compensation moments and performs R2s
%estimation. 
%  
% v1: Change to a global loop oover length(par.acuqistions)
%
%   Author: Martin Soellradl
%   Department of Neurology, Medical University of Graz, Graz, Austria
%   email:martin.soellradl@medunigraz.at
%   Website: http://www.neuroimaging.com
%   Januray 2020; Last revision: 20-Spetember-2020
   
    %% Convert the dcm files from to nifi (to get nifi header) and load 
    %headers 
   
    if exist([par.src_nii, 'dcm_header_short.mat'], 'file') ~= 2
        dicm2nii(par.src_dcm, par.src_nii); 
    end
    %load header
    tmp = load([par.src_nii, 'dcm_header_short.mat']);
    dcm_header = tmp.dcm_header_short;
    
    data = []; 

    %index of desired header 
    for i= 1:length(par.acquisitions); 
        

        %create a nifti template 
        tmp = load_untouch_nii([par.src_nii, 'Example_5_Gc_split_max_60deg.nii.gz']); 
        %store header
        tmp.hdr.dime.datatype = 16; 
        tmp.hdr.dime.bitpix = 32;   
        tmp.hdr.dime.dim(5) = tmp.hdr.dime.dim(5)-1; %reduce due to navigator
        nii_template_4D = tmp;

        nii_headers  = nii_template_4D; 

        tmp.hdr.dime.dim([1 5])=[4 1];
        nii_template = tmp; 


         %% Load and combine raw data with and without navigator echoes

        navi = {'navi_off', 'navi_on'};
        
        fwhm_str = num2str(dcm_header.SliceThickness);
        fwhm_str(fwhm_str == '.') = '_'; 
        fwhm_str = ['fwhm_', fwhm_str, '_mm']; 


        flip_angles = ['alpha_', num2str(dcm_header.FlipAngle), 'deg']; 

        S_uncorr_path = [par.path_pwd, '/results/R2s_zshim/', par.meas_id, '/', par.acquisitions{i}, '/', navi{1},'/', flip_angles, '/', fwhm_str, '/raw_images/']; 
        S_corr_path =  [par.path_pwd, '/results/R2s_zshim/', par.meas_id, '/', par.acquisitions{i}, '/', navi{2},'/', flip_angles, '/', fwhm_str, '/raw_images/']; 

        bExample = 1;
        [Scomb_corr, Scomb_raw  ] = prepareAndCorrectRawData( par.dat_path{i}, S_uncorr_path, S_corr_path, par.bdetrend, nii_headers, bExample);

        S_navi_on = double(Scomb_corr); 
        S_navi_off = double(Scomb_raw);

        %% Check if B1-map is avaiable 

        if par.bB1_map == 1
                tmp = load(par.B1_path);
                B1_map = imrotate(flip(tmp.flipAngleMap/100,2),-90);
                [Ny, Nx, Nz, Nte] = size(S_navi_on); 
                B1_map = imresize(B1_map, [Ny, Nx]); 
        end
    

        %% Start R2* pipeline 


        te = dcm_header.EchoTime; te = te(1:end-1); %delete navigator echo time
        z0 = dcm_header.SliceThickness;

        opts = []; 
        opts.B0_method  = 'phase_diff_B0';
        opts.bbet = par.bbet; 
        opts.sgl_sag_slc = par.sgl_sag_slc; 
        opts.bprelude = 1;
        opts.nii_template = nii_template; 
        opts.nii_template_4D = nii_template_4D;


        if isfield(par, 'z0_vary')
            z0_vary = par.z0_vary; 
        else
            z0_vary = dcm_header.SliceThickness; 
        end
        
        pulse =par.pulse_c{i};
        pulse.alpha =  dcm_header.FlipAngle; 

        %check if parameters only with navigator data should be estimated
        if par.navi_only == 1
          idx_navi = 2; 
        else
          idx_navi = [1,2]; 
        end



        %Calculate shim moments
        if strcmp( par.zShim.pattern{i}, '+A_-A_0_+B_-B')
             if isnan(par.zShim.Gc{i}) %Value is NaN if a slice specific pattern is selected
                 Gz_comp =  dlmread(par.zShim.Gz_comp_path);   
             else 
                 Gz_comp = par.zShim.Gc{i}.*ones(size(S_navi_on,3), 2)*1E-3; %mT/m
                 Gz_comp(:,2) =  -1*Gz_comp(:,2);
             end
        else
             Gz_comp =  dlmread(par.zShim.Gz_comp_path);    
        end
        
        %Based on par.zShim.pattern{i} the shim moment for before each is
        %returned. 
        [ MzShim, GzShim] = getGzShimPattern( Gz_comp, te, par.zShim.pattern{i}, par.zShim.start_echo); 


        for k= idx_navi

            if k==1
               mag = abs(S_navi_off); 
               phase = angle(S_navi_off); 
            else
               mag = abs(S_navi_on); 
               phase = angle(S_navi_on); 
            end

            for l=1:length(z0_vary)   

              fwhm_str = num2str(z0_vary(l));

              fwhm_str(fwhm_str == '.') = '_'; 
              fwhm_str = ['fwhm_', fwhm_str, '_mm']; 

        
              path_results = [par.path_pwd, '/results/R2s_zshim/', par.meas_id, '/', par.acquisitions{i}, '/', navi{k},'/', flip_angles, '/', fwhm_str]; 


              if exist(path_results, 'dir') == 0
                  mkdir(path_results); 
              end
              path_src_nii = [path_results, '/raw_images/mag_', navi{k}, '.nii.gz']

              file_id = [par.acquisitions{i}, '_', flip_angles];

               if par.bB1_map
                    results = pipelineR2sCorrectionWithPulseShapeAndZShimming_v1(mag, phase, te, MzShim, pulse, file_id, path_results, path_src_nii, z0_vary, opts, B1_map);
               else
                    results = pipelineR2sCorrectionWithPulseShapeAndZShimming_v1(mag, phase, te, MzShim, pulse, file_id, path_results, path_src_nii, z0_vary, opts);
               end
               
               %Estiamte R2* from a monoexponential for the refrence image
               %w/o z-shim
               if strcmp( par.acquisitions{i}, 'zShim_off')
            
                 [ res ] = performMonoexponentialFit(mag, results.mask, te,file_id, path_results, opts)
                 results.R2s_mono = res.R2s_mono;  
                 
               end

            end
            
                      
            

            %store results in original struct 
            data = setfield(data, par.acquisitions{i}, 'results', navi{k}, flip_angles,results); 
            data = setfield(data, par.acquisitions{i}, 'results', navi{k}, flip_angles,'dcm_header', dcm_header);
            data = setfield(data, par.acquisitions{i}, 'results', navi{k}, flip_angles,'mag', mag);
            data = setfield(data, par.acquisitions{i}, 'results', navi{k}, flip_angles,'phase', phase);
         
        end 
    end


end

