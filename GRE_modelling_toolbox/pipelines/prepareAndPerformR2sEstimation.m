function [ data ] = prepareAndPerformR2sEstimation( par )
%PREPAREANDPERFORMMWFESTIMATION Script loads SIEMENS raw data with
%navigator echo, performs a corrected coil combination and calls the
%pipeline for R2* esimation. 
%
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
    header_all = load([par.src_nii, 'dcmHeaders']);
    header_all = header_all.h; 

    %get file names from header 
    file_names = fieldnames(header_all);

    
    data = [];    
    
    %index of desired header 
    for i= 1:length(par.acquisitions); 

        %create a nifti template for saving estimated maps and corrected
        %raw data 
        tmp = load_untouch_nii([par.src_nii, file_names{par.idx_nii(i)}, '.nii.gz']); 
        
        %store header information
        tmp.hdr.dime.datatype = 16; 
        tmp.hdr.dime.bitpix = 32;   
        tmp.hdr.dime.dim(5) = tmp.hdr.dime.dim(5)-1; %reduce dimension due to the navigator echo
        nii_template_4D = tmp;

        %4D template
        nii_headers  = nii_template_4D; 

        %3D template
        tmp.hdr.dime.dim([1 5])=[4 1];
        nii_template = tmp; 


         %% Load and combine raw data with and without navigator echoes

        dcm_header = getfield(header_all, file_names{par.idx_nii(i)});

        navi = {'navi_off', 'navi_on'};

        fwhm_str = num2str(dcm_header.SliceThickness);
        fwhm_str(fwhm_str == '.') = '_'; 
        fwhm_str = ['fwhm_', fwhm_str, '_mm']; 


        flip_angles = ['alpha_', num2str(dcm_header.FlipAngle), 'deg']; 

        path_results = [par.path_pwd, '/results/R2s/', par.meas_id, '/', par.acquisitions{i} ,'/']; 
        S_uncorr_path = [par.path_pwd, '/results/R2s/', par.meas_id, '/', par.acquisitions{i}, '/', navi{1},'/', flip_angles, '/', fwhm_str, '/raw_images/']; 
        S_corr_path =  [par.path_pwd, '/results/R2s/', par.meas_id, '/', par.acquisitions{i}, '/', navi{2},'/', flip_angles, '/', fwhm_str, '/raw_images/']; 

        %Raw data is laoded and corrected with the navigator echo.
        %Corrected raw data is combined with the method proposed by Luo et
        %al
        [Scomb_corr, Scomb_raw  ] = prepareAndCorrectRawData( par.dat_path{i}, S_uncorr_path, S_corr_path, par.bdetrend, nii_headers);

        S_navi_on = double(Scomb_corr); 
        S_navi_off = double(Scomb_raw);

        %% Check if B1-map is avaiable 

        if par.bB1_map == 1

            if par.bB1_map_new==1
                tmp = load(par.B1_path);
                B1_map = imrotate(flip(tmp.flipAngleMap/100,2),-90);
                [Ny, Nx, Nz, Nte] = size(S_navi_on); 
                B1_map = imresize(B1_map, [Ny, Nx]); 
            else 
                tmp = load(par.B1_path);
                B1_map = flip(tmp.flipAngleMap/100);
                [Ny, Nx, Nz, Nte] = size(S_navi_on); 
                 B1_map = imresize(B1_map, [Ny, Nx]); 
            end
        end
        
        %% Start R2* pipeline 

        te = dcm_header.EchoTime; te = te(1:end-1); %delete navigator echo time
        z0 = dcm_header.SliceThickness;

        %set options for R2* esimation
        opts = []; 
        opts.B0_method  = 'phase_diff_B0';
        opts.bbet = par.bbet; 
        opts.sgl_sag_slc = par.sgl_sag_slc; 
        opts.bprelude = 1;
        opts.nii_template = nii_template; 
        opts.nii_template_4D = nii_template_4D;
        opts.Gphi0z = par.Gphi0z; 


        %regularization option. For in-vivo measurements the nominal slice
        %thickness leads to good results. 
        if isfield(par, 'z0_vary')
            z0_vary = par.z0_vary; 
        else
            z0_vary = dcm_header.SliceThickness; 
        end




       %extract pulse information
       pulse =par.pulse_c{i};
       pulse.alpha =  dcm_header.FlipAngle; 
       
       %check if parameters only with navigator data should be estimated
       if par.navi_only == 1
          idx_navi = 2; 
       else
          idx_navi = [1,2]; 
       end

       
       for k= idx_navi

           %set mag and phase data (either with or without navigator
           %correction)
           if k==1
               mag = abs(S_navi_off); 
               phase = angle(S_navi_off); 
           else
               mag = abs(S_navi_on); 
               phase = angle(S_navi_on); 
           end
        
           for l=1:length(z0_vary)

              %create string for file storage
              fwhm_str = []; 
              fwhm_str = num2str(z0_vary(l));
              fwhm_str(fwhm_str == '.') = '_'; 
              fwhm_str = ['fwhm_', fwhm_str, '_mm']; 

              if par.register_maps == 1 && i==1
                  path_results = [par.path_pwd, '/results/R2s/', par.meas_id, '/', par.acquisitions{i}, '_reg/', navi{k},'/', flip_angles, '/', fwhm_str]; 
              else
                  path_results = [par.path_pwd, '/results/R2s/', par.meas_id, '/', par.acquisitions{i}, '/', navi{k},'/', flip_angles, '/', fwhm_str]; 
              end

              if exist(path_results, 'dir') == 0
                  mkdir(path_results); 
              end
              path_src_nii = [path_results, '/raw_images/mag_', navi{k}, '.nii.gz']

              file_id = [par.acquisitions{i}, '_', flip_angles];

              z0_vary = z0; 

              %start R2s estimation depending on paramters
              if par.register_maps == 1 && i==1
                 file_id = [par.acquisitions{i}, '_reg_', flip_angles];
                 if par.bB1_map == 1
                     results =  pipelineR2sCorrectionWithPulseShape_v5_loadRegMaps(mag, phase, te, pulse, file_id, path_results, path_src_nii, z0, z0_vary,opts, B1_map);
                 else
                    results =  pipelineR2sCorrectionWithPulseShape_v5_loadRegMaps(mag, phase, te, pulse, file_id, path_results, path_src_nii, z0, z0_vary,opts);
                 end
              else
                 if par.bB1_map == 1
                    results =  pipelineR2sCorrectionWithPulseShape_v5(mag, phase, te, pulse, file_id, path_results, path_src_nii, z0, z0_vary(l),opts, B1_map);
                 else
                    results =  pipelineR2sCorrectionWithPulseShape_v5(mag, phase, te, pulse, file_id, path_results, path_src_nii, z0, z0_vary(l),opts);
                 end
             end
             
             
              %Segment WM from MP2RAGE image and register mask on mGRE imags 
              if par.bSienna == 1
                 fwhm_str = num2str(z0);
                 fwhm_str(fwhm_str == '.') = '_'; 
                 fwhm_str = ['fwhm_', fwhm_str, '_mm']; 
                 path_results = [par.path_pwd, '/results/R2s/', par.meas_id, '/', par.acquisitions{i}, '/', navi{k},'/', flip_angles, '/', fwhm_str]; 
                 par.path_results = path_results; 
                 par.file_id = file_id;
                 performSIENNAandRegisterToGRE( par) 
              end
   
             
          end
          

          %store results in original struct 
          data = setfield(data, par.acquisitions{i}, 'results', navi{k}, flip_angles,results); 
          data = setfield(data, par.acquisitions{i}, 'results', navi{k}, flip_angles,'dcm_header', dcm_header);
          data = setfield(data, par.acquisitions{i}, 'results', navi{k}, flip_angles,'mag', mag);
          if par.bB1_map == 1
            data = setfield(data, par.acquisitions{i}, 'results', navi{k}, flip_angles,'B1_map', B1_map);
          end
          data = setfield(data, par.acquisitions{i}, 'results', navi{k}, flip_angles,'phase', phase);
          data = setfield(data, par.acquisitions{i}, 'results', navi{k}, flip_angles,'nii_template', nii_template);
          data = setfield(data, par.acquisitions{i}, 'results', navi{k}, flip_angles,'nii_template_4D', nii_template_4D);
       end
    end

    % Register both acquistions 

    if par.register_maps == 1
        ref_path =  [par.src_nii, file_names{par.idx_nii(2)}]; 
        in_path =  [par.src_nii, file_names{par.idx_nii(1)}]; 
        %destiation folder for registered images 
        dst_folder_path =  [par.path_pwd, '/results/R2s/', par.meas_id, '/', par.acquisitions{1}];

        if ~exist([dst_folder_path, '_reg'], 'dir')
            registrationResults(ref_path, in_path, dst_folder_path);
        end
    end
end

