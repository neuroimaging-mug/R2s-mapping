function [ data ] = prepareAndPerformMWFEstimation_v2( par )
%PREPAREANDPERFORMMWFESTIMATION Summary of this function goes here
%   Detailed explanation goes here. 
    %v1: We add the regulariaztion paramter for the MWF estimation as
    %addtional paramter. 
    %v2: Online version with anonymized data. 
    
%   Author: Martin Soellradl
%   Department of Neurology, Medical University of Graz, Graz, Austria
%   email:martin.soellradl@medunigraz.at
%   Website: http://www.neuroimaging.com
%   Januray 2020; Last revision: 02-Januray-2020



    %% Convert the dcm files from to nifi (to get nifi header) and load 
    %headers 
       
    if exist([par.src_nii, 'dcmHeaders.mat'], 'file') ~= 2
        dicm2nii(par.src_dcm, par.src_nii); 
    end
    %load header
    load([par.src_nii, 'dcm_header_short']);

    %index of desired header 
    for i=1:length(par.idx_nii); 

        %create a nifti template 
        tmp = load_untouch_nii( par.nii_file ); 
        %store header
        tmp.hdr.dime.datatype = 16; 
        tmp.hdr.dime.bitpix = 32;   
        tmp.hdr.dime.dim(5) = tmp.hdr.dime.dim(5)-1; %reduce due to navigator
        nii_template_4D = tmp;

        nii_headers{i}  = nii_template_4D; 

    end
    tmp.hdr.dime.dim([1 5])=[4 1];
    nii_template = tmp; 
    
    
    %% Load raw data, perform navi correction and combine data 
    
    
    fwhm_str = ['fwhm_4_mm']; 
    flip_angles{1} = ['alpha_', num2str(dcm_header.FlipAngle), 'deg']; 
    navi = {'navi_off', 'navi_on'};
    dst_raw =  [par.path_pwd, '/results/MWF/', par.meas_id, '/',   par.acquisition , '/', navi{1},'/', flip_angles{1}, '/', fwhm_str, '/raw_images/']; 
    dst_navi = [par.path_pwd, '/results/MWF/', par.meas_id, '/',   par.acquisition , '/', navi{2},'/', flip_angles{1}, '/', fwhm_str, '/raw_images/']; 
 
    res_avg = averageNavigatorData( par.dat_avg_path, nii_headers, dst_raw, dst_navi,par.bdetrend, par.bloadExample);
    
    %% define pulses 

    %myRF fast
    pulse.type = 'sinc-hanning'; 
    pulse.Tpulse = 1;               %ms 
    pulse.BWT = 2;                %-
    pulse.k_pulse = 56.6132;         %�T
    polarity = {'positive'}; 

    
    %% Check if B1-map is avaiable 
    
    if par.bB1_map == 1
        
        if strcmp(par.meas_id, 'subject_6_Jo')
            tmp = load(par.B1_path);
            B1_map = flip(tmp.flipAngleMap/100);
            [Ny, Nx, Nz, Nte] = size(res_avg.mag_avg_navi); 
            B1_map = imresize(B1_map, [Ny, Nx]); 
            %replicate last slice to get the same number of slices as for the magnitdue
            %image
            B1_map(:,:,29) = B1_map(:,:,28); 
        elseif par.B1_new == 1 %we have changed acquistion. Other B1-maps were not so good
            [Ny, Nx, Nz, Nte] = size(res_avg.mag_avg_navi); 
            tmp = load(par.B1_path);
            B1_map = imrotate(flip(tmp.flipAngleMap/100,2),-90);
            B1_map = imresize(B1_map, [Ny, Nx]); 
        else
            tmp = load(par.B1_path);
            B1_map = imrotate(flip(tmp.flipAngleMap/100),-90);
            mag = res_avg.mag_avg_navi(:,:,:,1);     
            imagine(mag, B1_map); 
        end
    end
    
 

    %% Start R2* pipeline 

    opts = []; 
    opts.B0_method  = 'phase_diff_B0';
    opts.bbet = 1; 
    opts.bprelude = 0;
    opts.nii_template = nii_template; 
    opts.nii_template_4D = nii_template_4D;
    opts.AlonsoCorr = 0; 
    z0_vary = 4.0; %mm
    opts.beta_reg = par.beta_reg; 

    data = []; 

   for k=2 %just navi on data 

       if k==1
           mag = res_avg.mag_avg_raw; 
           phase = res_avg.phase_avg_raw; 
       else
           mag = res_avg.mag_avg_navi; 
           phase = res_avg.phase_avg_navi;
       end

       for j=1:length(flip_angles)

              path_results = [par.path_pwd, '/results/MWF/', par.meas_id, '/',  par.acquisition, '/', navi{k},'/', flip_angles{j}, '/', fwhm_str]; 
              if exist(path_results, 'dir') == 0
                  mkdir(path_results); 
              end
              path_src_nii = [path_results, '/raw_images/mag_', navi{k}, '.nii.gz']

              file_id = [par.acquisition, '_', flip_angles{j}];

              te = dcm_header.EchoTime; te = te(1:end-1); %delete navigator echo time
              z0 = dcm_header.SliceThickness;
              pulse.alpha =  dcm_header.FlipAngle; 
              pulse.GsPolarity = polarity{1}; 

              if par.bB1_map == 1
                results = pipelineFieldCorrectionWithPulseShapeForMWF_NNLS_MERA_v2(mag, phase, te, pulse, file_id, path_results, path_src_nii, z0, opts, B1_map);
              else
                 results = pipelineFieldCorrectionWithPulseShapeForMWF_NNLS_MERA_v2(mag, phase, te, pulse, file_id, path_results, path_src_nii, z0, opts);
              end
       
              %store results in original struct 
              data = setfield(data, par.acquisition, 'results', navi{k}, flip_angles{j},results); 
              data = setfield(data, par.acquisition, 'results', navi{k}, flip_angles{j},'dcm_header', dcm_header);
              data = setfield(data, par.acquisition, 'results', navi{k}, flip_angles{j},'mag', mag);
              data = setfield(data, par.acquisition, 'results', navi{k}, flip_angles{j},'phase', phase);
       end 
   end





end

