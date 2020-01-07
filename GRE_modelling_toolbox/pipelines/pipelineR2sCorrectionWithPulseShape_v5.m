function [ res ] = pipelineR2sCorrectionWithPulseShape_v5(mag, phase, te, pulse, file_id, path_results, path_src_nii, z0, z0_vary, opts, B1_map)
%PIPELINER2SCORRECTIONWITHPULSESHAPE Pipeline for R2s correction
% v1: Change of the gradient calcuation. Instead of using the linear
% regression approach with 26 neighbors fo a voxel the MATLAB gradient
% funciton is used. 
% v2: Change in the calcuation of Fn with the Bloch equations. A new
% function that includes several options (with, w/0 
% B1 / slice broading/thinning) is  included. 
% v5: The slice thickness of the foward simulation is changed to z0_vary to
% allow a change of the slice slecection gradient.If z0_vary = z0 calcuations 
%are identical to the previous calculatios. 
%
%   Author: Martin Soellradl
%   Department of Neurology, Medical University of Graz, Graz, Austria
%   email:martin.soellradl@medunigraz.at
%   Website: http://www.neuroimaging.com
%   Januray 2020; Last revision: 02-Januray-2020

    %Check if B1-map is available 
    bB1_map = 0; 
    if nargin == 11
       bB1_map = 1;  
    end

 %% store mag and phase files of the the echoes for further processing
    
    if opts.bbet == 1 || opts.bprelude == 1
            disp(['Storing mag and phase files for prelude and mask ', file_id, '...']); 
            path=[path_results, '/raw_images/'];
            if ~exist(path, 'dir')
                mkdir(path); 
            end
            
            mag_path =  [path, 'mag_', file_id, '.nii.gz'];
            phase_path = [path, 'phase_', file_id,'.nii.gz'];
            if ~exist(mag_path, 'file')
              
                nii = opts.nii_template_4D;
                nii.img = mag;
                save_untouch_nii(nii, mag_path); 
                nii.img = phase;
                save_untouch_nii(nii,phase_path); 

                for j=1:3 %size(mag{i},4); 
                    nii = opts.nii_template;
                    nii.img = mag(:,:,:,j);
                    save_untouch_nii(nii, [path, 'mag_echo', num2str(j), '_', file_id, '.nii.gz']); 
                    nii.img = phase(:,:,:,j);
                    save_untouch_nii(nii, [path, 'phase_echo', num2str(j), '_', file_id, '.nii.gz']); 
                end     
            end
    end
    
    
    
    %% Create Mask, either with BET or by simple thresholding 
    
    if opts.bbet == 1
        
         path = [path_results, '/Mask/'];
         if ~exist(path, 'dir')
            mkdir(path); 
         end
         gre_bet_path =  [path, 'bet_mag_echo1_', file_id, '.nii'];
         if exist(gre_bet_path, 'file') == 0
              %bet the gradient echo imagse (flirt works best for bet images)
              gre_path = [path_results '/raw_images/', 'mag_echo1_', file_id, '.nii.gz'];

              opts_bet = '-m -f 0.15 -R -S -B -Z';
              bet_call = ['bet ', gre_path, ' ', gre_bet_path, ' ', opts_bet]; 
              system(bet_call); 
         end
         %load bet-mask 
         mask_path = [path, 'bet_mag_echo1_', file_id, '_mask.nii.gz']; 
         tmp = load_untouch_nii(mask_path); 
         mask = double(tmp.img); 
         
    else % Create simple mask 
            path=[path_results, '/Mask/']; 
            if exist(path, 'dir') == 0
                mkdir(path);   
            end
            mask_path = [path, 'Mask_', file_id, '.nii']; 
            if exist(mask_path, 'file') == 0
                %create a simple mask 
                mask_tmp =  mag(:,:,:,1);
                th = 0.025*max(mask_tmp(:));
                mask_tmp(mask_tmp < th) = 0;  
                mask_tmp(mask_tmp >= th) = 1; 
                mask = mask_tmp; 

                %save mask 
                nii = opts.nii_template;
                nii.img = mask; 
                save_untouch_nii(nii, [path, 'Mask_', file_id, '.nii']); 
            else
                tmp = load_untouch_nii(mask_path); 
                mask = double(tmp.img);     
            end     
    end
    res.mask = mask; 
    
    
   
    
    %% If prelude flag is set data is unwrapped else it is is assumed that
    % the phase is already unrwapped!
    
    if opts.bprelude == 1
        %prelude -p <phasevol> -a <absvol> -o <unwrappedphase> [options]
          disp(['Unwrapping phase with prelude ', file_id, '...']); 
          
          path=[path_results, '/B0_map/prelude/unwrapped/'];
          if ~exist(path, 'dir')
              mkdir(path); 
          end

          unwrapped = [path, 'unwrapped_phase_', file_id, '.nii.gz']; 
          if ~exist(unwrapped, 'file')
              %path to mag and phase (4D)
              path = [path_results, '/raw_images/']; 
              phasevol = [path, 'phase_', file_id, '.nii.gz']; 
              absvol = [path,  'mag_', file_id, '.nii.gz'];

              %destination unwrapped phase 
              path = [path_results, '/B0_map/prelude/unwrapped/']; 
              mkdir(path);


              prelude_call = ['prelude ', '-p ', phasevol, ' -a ', absvol, ...
                  ' -o ',  unwrapped, ' -s']; 
              
              system(prelude_call); 

              %load unwrapped mask 
              tmp = load_untouch_nii(unwrapped); 
              phase_prelude = double(tmp.img); 
          else 
              tmp = load_untouch_nii(unwrapped); 
              phase_prelude = double(tmp.img); 
          end
          res.phase_prelude = phase_prelude;       
    end



  %% Create B0-Map 

  
 
   switch opts.B0_method
       case 'lin_fit_B0'
            disp( ['Linear  fit of B0-map ', file_id, '...']); 
            path=[path_results, '/B0_map/B0_lin_fit/'];
            if ~exist(path, 'dir')
                mkdir(path); 
            end
            phi0_path = [path, 'phi0_linfit_', file_id, '.nii']; 
            dw0_path = [path, 'dw0_linfit_',  file_id, '.nii'];
            
            if ~exist(phi0_path, 'file'); 
                te_selected = te(1:4); 

                bUnwrap = 0; 
                [ dw0, phi0] = LinearFitOfB0map(phase_prelude,mask, te_selected, 0, bUnwrap);

                nii = opts.nii_template;
                nii.img = phi0; 
                save_untouch_nii(nii, phi0_path); 
                nii.img = dw0; 
                save_untouch_nii(nii, dw0_path); 
            else
               tmp = load_untouch_nii(phi0_path); 
               phi0 = double(tmp.img); 
               tmp = load_untouch_nii(dw0_path); 
               dw0 = double(tmp.img); 
            end
   
            res.dw0 = dw0; 
            res.ph0 = phi0; 
       case 'cplx_fit_B0'
            disp( ['Complex fit of B0-map ', file_id, '...']); 
            
            path=[path_results, '/B0_map/B0_cplx_fit/'];
            if ~exist(path, 'dir')
                mkdir(path); 
            end
            
            phi0_path = [path, 'phi0_cplx_', file_id, '.nii']; 
            dw0_path = [path, 'dw0_cplx_',  file_id, '.nii'];
            if ~exist(phi0_path, 'file'); 
                te_selected = te(1:10); 

                [ dw0, phi0] = CplxFitOfB0map(phase,mask, te_selected, 0);
                
                nii = opts.nii_template;
                nii.img = phi0; 
                save_untouch_nii(nii, phi0_path); 
                nii.img = dw0; 
                save_untouch_nii(nii, dw0_path); 
            else
               tmp = load_untouch_nii(nii, phi0_path); 
               phi0 = double(tmp.img); 
               tmp = load_untouch_nii(nii, dw0_path); 
               dw0 = double(tmp.img);  
                
            end
            res.dw0 = dw0; 
            res.ph0 = phi0;   
       case 'phase_diff_B0'
    
            phi0 = zeros(size(mag,1), size(mag,2), size(mag,3)); 
            path=[path_results, '/B0_map/B0_from_prelude/Echo12/']; 
           
            if ~exist(path, 'dir')
                mkdir(path); 
            end  
            
            path_B0_12 = [path, 'B0_from_prelude_12_', file_id, '.nii'];
            path_phi0_12 = [path, 'phi0_from_prelude_12_', file_id, '.nii']
            if ~exist(path_B0_12, 'file') || ~exist(path_phi0_12, 'file')
          
                dphi =  phase_prelude(:,:,:,2) - phase_prelude(:,:,:,1);
                dTE = (te(2) - te(1))*1E-3; 
                dw0_12 = dphi./dTE; 
                phi0_12 = phase_prelude(:,:,:,2) - dw0_12*te(2)*1E-3;

                nii = opts.nii_template; 
                nii.img =  dw0_12; 
                save_untouch_nii(nii, path_B0_12);
                nii.img =  phi0_12; 
                save_untouch_nii(nii, path_phi0_12);
            else
                tmp = load_untouch_nii([path_B0_12]); 
                dw0_12 = double(tmp.img); 
                tmp = load_untouch_nii(path_phi0_12); 
                phi0_12 = double(tmp.img); 
            end
            
            path=[path_results, '/B0_map/B0_from_prelude/Echo13/'];
            if ~exist(path, 'dir')
                mkdir(path); 
            end
            
            path_B0_13 = [path, 'B0_from_prelude_13_', file_id, '.nii'];
            path_phi0_13 = [path, 'phi0_from_prelude_13_', file_id, '.nii'];
            if ~exist(path_B0_13, 'file') ||  ~exist(path_phi0_13, 'file')
                dTE =  (te(3) -  te(1))*1E-3;
                dphi = phase_prelude(:,:,:,3) - phase_prelude(:,:,:,1);
                dw0_13 = dphi./dTE; 
                phi0_13 = phase_prelude(:,:,:,2) - dw0_13*te(2)*1E-3;
                
                
                nii = opts.nii_template; 
                nii.img =  dw0_13; 
                save_untouch_nii(nii, path_B0_13);
                nii.img =  phi0_13; 
                save_untouch_nii(nii, path_phi0_13);
            else
                tmp = load_untouch_nii(path_B0_13); 
                dw0_13 = double(tmp.img); 
                tmp = load_untouch_nii(path_phi0_13); 
                phi0_13 = double(tmp.img); 
            end
            
            dw0 = dw0_12; 
            phi0 = phi0_12; 
            res.dw0_12 = dw0_12; 
            res.dw0_13 = dw0_13;
            res.phi0_12 = phi0_12;
            res.phi0_13 = phi0_13;
   end
   
    
    %% Estimate gradients with finite differencs 
    
    
    disp( ['Gradient estimation B0-map ', file_id,  '...']); 
    path=[path_results '/Gradient_maps/', opts.B0_method, '/']; 
    if ~exist(path, 'dir')
            mkdir(path); 
    end
    
    Gx_path = [path, 'Gx_finite_', opts.B0_method,'_', file_id , '.nii'];
    Gy_path = [path, 'Gy_finite_', opts.B0_method,'_', file_id , '.nii'];
    Gz_path = [path, 'Gz_finite_', opts.B0_method,'_', file_id , '.nii'];
    Gmag_path = [path, 'Gmag_finite_', opts.B0_method,'_', file_id , '.nii'];
    
    if ~exist(Gz_path, 'file')
        
        [ Gx, Gy, Gz] = gradient(dw0);
        Gx = Gx.*mask; 
        Gy = Gy.*mask; 
        Gz = Gz.*mask; 
        
        Gmag = sqrt(Gx.^2 + Gy.^2 + Gz.^2); 

        mkdir(path); 
        nii = opts.nii_template; 

        nii.img = Gx; 
        save_untouch_nii(nii, Gx_path);
        nii.img = Gy; 
        save_untouch_nii(nii, Gy_path);
        nii.img = Gz; 
        save_untouch_nii(nii, Gz_path);
        nii.img =  Gmag; 
        save_untouch_nii(nii, Gmag_path);

    else
        tmp = load_untouch_nii(Gx_path); 
        Gx = double(tmp.img); 
        tmp = load_untouch_nii(Gy_path); 
        Gy = double(tmp.img); 
        tmp = load_untouch_nii(Gz_path); 
        Gz = double(tmp.img); 
    end

    res.Gx = Gx;   
    res.Gy = Gy;   
    res.Gz = Gz;
    


    % ----------------------------------------------------------------
    % Perform a monoexponential fit of the data (model S1)
    % ----------------------------------------------------------------
            

    disp( ['Conventional monoexponential fit (S1) of the data ', file_id,  '...']); 
    path=[path_results '/R2s_monoexp/']; 
    if ~exist(path, 'dir')
            mkdir(path); 
    end

    R2s_corr_path = [path, 'R2s_mono_', file_id, '.nii'];
    A_corr_path = [path, 'A_mono_', file_id, '.nii'];
    res_norm_path = [path, 'resnorm_mono_', file_id, '.nii'];
    residuals_path = [path, 'residuals_mono_', file_id, '.nii'];

    if ~exist(R2s_corr_path, 'file')

        F_ones = ones(size(mag)); 
        [R2s_mono, A_mono, resnorm_mono, residuals_mono] ...
            = CorrectedR2sMapEstimationFnNonLinFit(mag, F_ones, te, mask);

        mkdir(path); 
        nii = opts.nii_template; 

        nii.img = R2s_mono; 
        save_untouch_nii(nii, R2s_corr_path);
        nii.img =  A_mono; 
        save_untouch_nii(nii, A_corr_path);
        nii.img =  resnorm_mono; 
        save_untouch_nii(nii, res_norm_path);
        nii.img =  residuals_mono; 
        save_untouch_nii(nii, residuals_path);

    else
        tmp = load_untouch_nii(R2s_corr_path); 
        R2s_mono = double(tmp.img); 
        tmp = load_untouch_nii(A_corr_path); 
        A_mono = double(tmp.img); 
        tmp = load_untouch_nii(res_norm_path); 
        resnorm_mono = double(tmp.img); 

    end
    res.R2s_mono = R2s_mono; 
    
    
   
    
%% Estimate Fn function for the given data with the Preibisch model 
   
    disp( ['Estimage Fn ', file_id,  '...']);
    
    gamma = 267.51;                      % gyromagnetic ratio 42.57*2*pi in MHz/T
                             
    path_wSlcCorr=[path_results '/Fn_maps/wSlcCorr/']; 
    path_preibisch=[path_results '/Fn_maps/Preibisch_approach/']; 
    if ~exist(path_wSlcCorr, 'dir')
            mkdir(path_wSlcCorr); 
    end
    if ~exist(path_preibisch, 'dir')
            mkdir(path_preibisch); 
    end
    
    Fn_wSlcCorr_path = [path_wSlcCorr, 'Fn_wSlcCorr_',  file_id , '.nii'];
    Fn_preibisch_path = [path_preibisch, 'Fn_preibisch_',  file_id , '.nii'];

    if ~exist(Fn_preibisch_path, 'file')
        
        Gsus = Gz/(gamma*z0); %convert from rad/(s*slice_thicness) to mT/m
        [ Fn_wSlcCorr, Fn_preibisch ] = estimateFnMapWtihRFShape( pulse, Gsus, te, z0_vary );

        nii = opts.nii_template_4D; 

        nii.img = Fn_wSlcCorr; 
        save_untouch_nii(nii, Fn_wSlcCorr_path);
        nii.img = Fn_preibisch; 
        save_untouch_nii(nii, Fn_preibisch_path);

    else
        tmp = load_untouch_nii(Fn_wSlcCorr_path); 
        Fn_wSlcCorr = double(tmp.img); 
        tmp = load_untouch_nii(Fn_preibisch_path); 
        Fn_preibisch = double(tmp.img); 

    end

    res.Fn_preibisch = Fn_preibisch;   
    res.Fn_wSlcCorr =  Fn_wSlcCorr;   
    
    

    % In case phi0 was also estimated from teh data (requires different
    % coil combination) the gradient Gphi0z can also be included. 
    if isfield(opts, 'Gphi0z');      
        if  opts.Gphi0z == 1           

            %% Estimate gradients from of phi0 if linear fit was performed
     
            disp( ['Gradient estimation from phi0 map ', file_id,  '...']); 
            path=[path_results '/Gradient_maps/', opts.B0_method, '/']; 
            if ~exist(path, 'dir')
                    mkdir(path); 
            end
        
            Gx_path = [path, 'G_phi0x_finite_', opts.B0_method,'_', file_id , '.nii'];
            Gy_path = [path, 'G_phi0y_finite_', opts.B0_method,'_', file_id , '.nii'];
            Gz_path = [path, 'G_phi0z_finite_', opts.B0_method,'_', file_id , '.nii'];
            Gmag_path = [path, 'Gmag_phi0_finite_', opts.B0_method,'_', file_id , '.nii'];
        
            if ~exist(Gz_path, 'file')
        
                [ G_phi0x, G_phi0y, G_phi0z] = gradient(phi0);
                
                G_phi0z = imgaussfilt3(G_phi0z,1); %filter gradient maps
                G_phi0x = imgaussfilt3(G_phi0x,1);
                G_phi0y = imgaussfilt3(G_phi0y,1);
                
                G_phi0x = G_phi0x.*mask; 
                G_phi0y = G_phi0y.*mask; 
                G_phi0z = G_phi0z.*mask; 
        
                Gmag = sqrt(G_phi0x.^2 + G_phi0y.^2 + G_phi0z.^2); 
        
                mkdir(path); 
                nii = opts.nii_template; 
        
                nii.img = G_phi0x; 
                save_untouch_nii(nii, Gx_path);
                nii.img = G_phi0y; 
                save_untouch_nii(nii, Gy_path);
                nii.img = G_phi0z; 
                save_untouch_nii(nii, Gz_path);
                nii.img =  Gmag; 
                save_untouch_nii(nii, Gmag_path);
        
            else
                tmp = load_untouch_nii(Gx_path); 
                G_phi0x = double(tmp.img); 
                tmp = load_untouch_nii(Gy_path); 
                G_phi0y = double(tmp.img); 
                tmp = load_untouch_nii(Gz_path); 
                G_phi0z = double(tmp.img); 
            end
        
            res.G_phi0x = G_phi0x;   
            res.G_phi0y = G_phi0y;   
            res.G_phi0z = G_phi0z;
            
    

            % ----------------------------------------------------------------
            % Estiamte Fn inclduing the gradient of phi0 without B1/lambda  in 
            % in z-direction 
            % ----------------------------------------------------------------
            
            
            disp( ['Estimation of Fn by solving the Bloch equations for the slice profile with phi0_z without B1 and slice correction ', file_id,  '...']); 
            path=[path_results '/Fn_maps/BlochNoB1_NoSlcCorr_Gphi0z/']; 
            if ~exist(path, 'dir')
                    mkdir(path); 
            end

            Fn_map_bloch_path = [path, 'Fn_map_BlochNoB1_noSlcCorr_Gphi0z', opts.B0_method,'_', file_id , '.nii'];
            if ~exist(Fn_map_bloch_path, 'file')

                bSlcCorr = 0; %slice correction due to the field gradient
                Gsus = Gz/(gamma*z0); %convert from rad/(s*slice_thicness) to mT/m
                G_phi0z_normalized = G_phi0z/z0; % rad/mm
                [ Fn_BlochNoB1NoSlcCorr_Gphi0z] =  estimateFnMapFromCplxSliceProfileWithBlochEquations_v1_withPhi0(pulse, Gsus, G_phi0z_normalized, z0_vary,  mask, te, bSlcCorr);

                nii = opts.nii_template_4D;  
                nii.img = Fn_BlochNoB1NoSlcCorr_Gphi0z; 
                save_untouch_nii(nii, Fn_map_bloch_path);
            else
                tmp = load_untouch_nii(Fn_map_bloch_path); 
                Fn_BlochNoB1NoSlcCorr_Gphi0z = double(tmp.img); 

            end
            res.Fn_BlochNoB1NoSlcCorr_Gphi0z = Fn_BlochNoB1NoSlcCorr_Gphi0z;


            %Estimate corrected R2s from Bloch equation with phase of pulse, without B1
            %and no slice correction. 
            disp( ['Corrected R2s estimation with Fn estimated from Bloch equations with phi0_z, no B1, no slice correction ', file_id,  '...']); 
            path=[path_results '/R2s_corr/BlochNoB1_NoSlcCorr_Gphi0z/']; 
            if ~exist(path, 'dir')
                    mkdir(path); 
            end

            R2s_corr_path = [path, 'R2s_BlochNoB1_noSlcCorr_Gphi0z', file_id, '.nii'];
            A_corr_path = [path, 'A_BlochNoB1_noSlcCorr_Gphi0z', file_id, '.nii'];
            res_norm_path = [path, 'resnorm_BlochNoB1_noSlcCorr_Gphi0z', file_id, '.nii'];
            residuals_path = [path, 'residuals_BlochNoB1_noSlcCorr_Gphi0z', file_id, '.nii'];

            %single slice option for faster calculation
            if opts.sgl_sag_slc == 1; 

               mask_opt = zeros(size(mag,1), size(mag,2), size(mag,3));
               mask_opt(:,size(mask,2)/2,:) = 1; 
               mask_opt = mask_opt.*double(mask); 
            else
               mask_opt = mask; 
            end

            if ~exist(R2s_corr_path, 'file')

                [R2s_corr_BlochNoB1NoSlcCorr_Gphi0z, A_corr_BlochNoB1NoSlcCorr_Gphi0z, resnorm_map_BlochNoB1NoSlcCorr_Gphi0z, residuals_map_BlochNoB1NoSlcCorr_Gphi0z] ...
                    = CorrectedR2sMapEstimationFnNonLinFit(mag, Fn_BlochNoB1NoSlcCorr_Gphi0z, te, mask_opt);

                mkdir(path); 
                nii = opts.nii_template; 

                nii.img = R2s_corr_BlochNoB1NoSlcCorr_Gphi0z; 
                save_untouch_nii(nii, R2s_corr_path);
                nii.img =  A_corr_BlochNoB1NoSlcCorr_Gphi0z; 
                save_untouch_nii(nii, A_corr_path);
                nii.img =  resnorm_map_BlochNoB1NoSlcCorr_Gphi0z; 
                save_untouch_nii(nii, res_norm_path);
                nii.img =  residuals_map_BlochNoB1NoSlcCorr_Gphi0z; 
                save_untouch_nii(nii, residuals_path);

            else
                tmp = load_untouch_nii(R2s_corr_path); 
                R2s_corr_BlochNoB1NoSlcCorr_Gphi0z = double(tmp.img); 
        %         tmp = load_untouch_nii(A_corr_path); 
        %         A_corr_BlochNoB1NoSlcCorr_Gphi0z = double(tmp.img); 
             %   tmp = load_untouch_nii(res_norm_path); 
              %  resnorm_map_BlochNoB1NoSlcCorr = double(tmp.img); 
            end

            res.R2s_corr_BlochNoB1NoSlcCorr_Gphi0z = R2s_corr_BlochNoB1NoSlcCorr_Gphi0z;   

        end
        
        
        
        if bB1_map == 1 && opts.Gphi0z
            

            % ----------------------------------------------------------------
            % Gphi0z With B1 and with slice correction
            % ----------------------------------------------------------------



            disp( ['Estimation of Fn by solving the Bloch equations for the slice profile including Gphi0z,  B1 errors and slice correction', file_id,  '...']); 
            path=[path_results '/Fn_maps/BlochWithB1_withSlcCorr_Gphi0z/']; 
            if ~exist(path, 'dir')
                    mkdir(path); 
            end

            Fn_map_bloch_path = [path, 'Fn_map_BlochWithB1_withSlcCorr_Gphi0z', opts.B0_method,'_', file_id , '.nii'];
            if ~exist(Fn_map_bloch_path, 'file')
                G_phi0z_normalized = G_phi0z/z0; % rad/mm
                bSlcCorr = 1; %slice correction due to the field gradient
                Gsus = Gz/(gamma*z0); %convert from rad/(s*slice_thicness) to mT/m
                [ Fn_BlochWithB1WithSlcCorr_Gphi0z] = ...
                    estimateFnMapFromCplxSliceProfileWithBlochEquations_v1_withPhi0(pulse, Gsus, G_phi0z_normalized, z0_vary,  mask, te, bSlcCorr, B1_map);


                nii = opts.nii_template_4D;  
                nii.img = Fn_BlochWithB1WithSlcCorr_Gphi0z; 
                save_untouch_nii(nii, Fn_map_bloch_path);
            else
                tmp = load_untouch_nii(Fn_map_bloch_path); 
                Fn_BlochWithB1WithSlcCorr_Gphi0z = double(tmp.img); 

            end
            res.Fn_BlochWithB1WithSlcCorr_Gphi0z = Fn_BlochWithB1WithSlcCorr_Gphi0z;



            %Estimate corrected R2s Fn esitmate from the phase including the
            % slice correction due to the field gradient of the pulse after excitation 
            disp( ['Corrected R2s estimation with Fn estimated from Bloch equations including Gphi0z, B1, slice correction, and phase of pulse', file_id,  '...']); 
            path=[path_results '/R2s_corr/BlochWithB1_withSlcCorr_Gphi0z/']; 
            if ~exist(path, 'dir')
                    mkdir(path); 
            end

            R2s_corr_path = [path, 'R2s_corr_BlochWithB1_withSlcCorr_Gphi0z', file_id, '.nii'];
            A_corr_path = [path, 'A_BlochWithB1_withSlcCorr_Gphi0z', file_id, '.nii'];
            res_norm_path = [path, 'resnorm_BlochWithB1_withSlcCorr_Gphi0z', file_id, '.nii'];
            residuals_path = [path, 'residuals_BlochWithB1_withSlcCorr_Gphi0z', file_id, '.nii'];

            if ~exist(R2s_corr_path, 'file')

                [R2s_corr_BlochB1_wSlcCorr_Gphi0z, A_corr_BlochB1_wSlcCorr_Gphi0z, resnorm_map_corr_BlochB1_wSlcCorr_Gphi0z, residuals_map_BlochB1_wSlcCorr_Gphi0z] ...
                    = CorrectedR2sMapEstimationFnNonLinFit(mag, Fn_BlochWithB1WithSlcCorr_Gphi0z, te, mask);

                mkdir(path); 
                nii = opts.nii_template; 

                nii.img = R2s_corr_BlochB1_wSlcCorr_Gphi0z; 
                save_untouch_nii(nii, R2s_corr_path);
                nii.img =  A_corr_BlochB1_wSlcCorr_Gphi0z; 
                save_untouch_nii(nii, A_corr_path);
                nii.img =  resnorm_map_corr_BlochB1_wSlcCorr_Gphi0z; 
                save_untouch_nii(nii, res_norm_path);
                nii.img =  residuals_map_BlochB1_wSlcCorr_Gphi0z; 
                save_untouch_nii(nii, residuals_path);

            else
                tmp = load_untouch_nii(R2s_corr_path); 
                R2s_corr_BlochB1_wSlcCorr_Gphi0z = double(tmp.img); 
%                 tmp = load_untouch_nii(A_corr_path); 
%                 A_corr_BlochB1_wSlcCorr_Gphi0z = double(tmp.img); 
%                 tmp = load_untouch_nii(res_norm_path); 
%                 resnorm_map_corr_BlochB1_wSlcCorr_Gphi0z = double(tmp.img); 
            end

            res.R2s_corr_BlochB1_wSlcCorr_Gphi0z = R2s_corr_BlochB1_wSlcCorr_Gphi0z; 

        end
             
            
            
            
            
    end


% 
%     ----------------------------------------------------------------
%     Estimate Fn with Bloch equations and only the mag of the slice 
%     profile (S2)
%     ----------------------------------------------------------------


    disp( ['Estimation of Fn by solving the Bloch equations for the magnitude slice profile ', file_id,  '...']); 
    path=[path_results '/Fn_maps/BlochMagOnly/']; 
    if ~exist(path, 'dir')
            mkdir(path); 
    end

    Fn_map_bloch_path = [path, 'Fn_map_BlochMagOnly_', opts.B0_method,'_', file_id , '.nii'];

    if ~exist(Fn_map_bloch_path, 'file')

        Gsus = Gz/(gamma*z0); %convert from rad/(s*slice_thicness) to mT/m
        [ Fn_BlochMagOnly] = estimateFnMapFromMagSliceProfileWithBlochEqn(pulse, Gsus, z0_vary,  mask, te);

        nii = opts.nii_template_4D;  
        nii.img = Fn_BlochMagOnly; 
        save_untouch_nii(nii, Fn_map_bloch_path);
    else
        tmp = load_untouch_nii(Fn_map_bloch_path); 
        Fn_BlochMagOnly = double(tmp.img); 

    end
    res.Fn_BlochMagOnly = Fn_BlochMagOnly;



    disp( ['Corrected R2s estimation with Fn estimated from Bloch equations by just using the magnitude (S2) ', file_id,  '...']); 
    path=[path_results '/R2s_corr/BlochMagOnly/']; 
    if ~exist(path, 'dir')
            mkdir(path); 
    end

    R2s_corr_path = [path, 'R2s_BlochMagOnly_', file_id, '.nii'];
    A_corr_path = [path, 'A_BlochMagOnly_', file_id, '.nii'];
    res_norm_path = [path, 'resnorm_BlochMagOnly_', file_id, '.nii'];
    residuals_path = [path, 'residuals_BlochMagOnly_', file_id, '.nii'];

    if ~exist(R2s_corr_path, 'file')

        [R2s_BlochMagOnly, A_BlochMagOnly, resnorm_BlochMagOnly, residuals_BlochMagOnly] ...
            = CorrectedR2sMapEstimationFnNonLinFit(mag, Fn_BlochMagOnly, te, mask);

        mkdir(path); 
        nii = opts.nii_template; 

        nii.img = R2s_BlochMagOnly; 
        save_untouch_nii(nii, R2s_corr_path);
        nii.img =  A_BlochMagOnly; 
        save_untouch_nii(nii, A_corr_path);
        nii.img =  resnorm_BlochMagOnly; 
        save_untouch_nii(nii, res_norm_path);
        nii.img =  residuals_BlochMagOnly; 
        save_untouch_nii(nii, residuals_path);

    else
        tmp = load_untouch_nii(R2s_corr_path); 
        R2s_BlochMagOnly = double(tmp.img); 
        tmp = load_untouch_nii(A_corr_path); 
            A_BlochMagOnly = double(tmp.img); 
            tmp = load_untouch_nii(res_norm_path); 
            resnorm_BlochMagOnly = double(tmp.img); 

    end
    res.R2s_BlochMagOnly = R2s_BlochMagOnly; 
    

    %% Estimate  Estimate Fn numerical by solving the Bloch equations 
    % including B1-errors
    if bB1_map == 1

      
        % ----------------------------------------------------------------
        % With B1 and with slice correction (S4)
        % ----------------------------------------------------------------
       

        disp( ['Estimation of Fn by solving the Bloch equations for the slice profile including B1 errors and slice correction', file_id,  '...']); 
        path=[path_results '/Fn_maps/BlochWithB1_withSlcCorr/']; 
        if ~exist(path, 'dir')
                mkdir(path); 
        end

        Fn_map_bloch_path = [path, 'Fn_map_BlochWithB1_withSlcCorr_', opts.B0_method,'_', file_id , '.nii'];
        if ~exist(Fn_map_bloch_path, 'file')

            bSlcCorr = 1; %slice correction due to the field gradient
            Gsus = Gz/(gamma*z0); %convert from rad/(s*slice_thicness) to mT/m
            [ Fn_BlochWithB1WithSlcCorr] =  estimateFnMapFromCplxSliceProfileWithBlochEquations_v1(pulse, Gsus, z0_vary,  mask, te, bSlcCorr, B1_map);
            

            nii = opts.nii_template_4D;  
            nii.img = Fn_BlochWithB1WithSlcCorr; 
            save_untouch_nii(nii, Fn_map_bloch_path);
        else
            tmp = load_untouch_nii(Fn_map_bloch_path); 
            Fn_BlochWithB1WithSlcCorr = double(tmp.img); 

        end
        res.Fn_BlochWithB1WithSlcCorr = Fn_BlochWithB1WithSlcCorr;
        
        
        
        %Estimate corrected R2s Fn esitmate from the phase including the
        % slice correction due to the field gradient of the pulse after excitation 
        disp( ['Corrected R2s estimation with Fn estimated from Bloch equations including B1, slice correction, and phase of pulse', file_id,  '...']); 
        path=[path_results '/R2s_corr/BlochWithB1_withSlcCorr/']; 
        if ~exist(path, 'dir')
                mkdir(path); 
        end

        R2s_corr_path = [path, 'R2s_corr_BlochWithB1_withSlcCorr_', file_id, '.nii'];
        A_corr_path = [path, 'A_BlochWithB1_withSlcCorr_', file_id, '.nii'];
        res_norm_path = [path, 'resnorm_BlochWithB1_withSlcCorr_', file_id, '.nii'];
        residuals_path = [path, 'residuals_BlochWithB1_withSlcCorr_', file_id, '.nii'];

        if ~exist(R2s_corr_path, 'file')

            [R2s_corr_BlochB1_wSlcCorr, A_corr_BlochB1_wSlcCorr, resnorm_map_corr_BlochB1_wSlcCorr, residuals_map_BlochB1_wSlcCorr] ...
                = CorrectedR2sMapEstimationFnNonLinFit(mag, Fn_BlochWithB1WithSlcCorr, te, mask);

            mkdir(path); 
            nii = opts.nii_template; 

            nii.img = R2s_corr_BlochB1_wSlcCorr; 
            save_untouch_nii(nii, R2s_corr_path);
            nii.img =  A_corr_BlochB1_wSlcCorr; 
            save_untouch_nii(nii, A_corr_path);
            nii.img =  resnorm_map_corr_BlochB1_wSlcCorr; 
            save_untouch_nii(nii, res_norm_path);
            nii.img =  residuals_map_BlochB1_wSlcCorr; 
            save_untouch_nii(nii, residuals_path);

        else
            tmp = load_untouch_nii(R2s_corr_path); 
            R2s_corr_BlochB1_wSlcCorr = double(tmp.img); 
            tmp = load_untouch_nii(A_corr_path); 
            A_corr_BlochB1_wSlcCorr = double(tmp.img); 
            tmp = load_untouch_nii(res_norm_path); 
            resnorm_map_corr_BlochB1_wSlcCorr = double(tmp.img); 
        end

        res.R2s_corr_BlochB1_wSlcCorr = R2s_corr_BlochB1_wSlcCorr; 
        %res.A_corr_BlochB1_wSlcCorr = A_corr_BlochB1_wSlcCorr; 
        %res.resnorm_map_corr_BlochB1_wSlcCorr = resnorm_map_corr_BlochB1_wSlcCorr; 
        

    end
        % ----------------------------------------------------------------
        % Without B1 and without slice correction (S3)
        % ----------------------------------------------------------------
       
          
        
        disp( ['Estimation of Fn by solving the Bloch equations for the slice profile without B1 and slice correction ', file_id,  '...']); 
        path=[path_results '/Fn_maps/BlochNoB1_NoSlcCorr/']; 
        if ~exist(path, 'dir')
                mkdir(path); 
        end

        Fn_map_bloch_path = [path, 'Fn_map_BlochNoB1_noSlcCorr_', opts.B0_method,'_', file_id , '.nii'];
        if ~exist(Fn_map_bloch_path, 'file')
            
            bSlcCorr = 0; %slice correction due to the field gradient
            Gsus = Gz/(gamma*z0); %convert from rad/(s*slice_thicness) to mT/m
            [ Fn_BlochNoB1NoSlcCorr] =  estimateFnMapFromCplxSliceProfileWithBlochEquations_v1(pulse, Gsus, z0_vary,  mask, te, bSlcCorr);
            
            nii = opts.nii_template_4D;  
            nii.img = Fn_BlochNoB1NoSlcCorr; 
            save_untouch_nii(nii, Fn_map_bloch_path);
        else
            tmp = load_untouch_nii(Fn_map_bloch_path); 
            Fn_BlochNoB1NoSlcCorr = double(tmp.img); 

        end
        res.Fn_BlochNoB1NoSlcCorr = Fn_BlochNoB1NoSlcCorr;
        
        
        %Estimate corrected R2s from Bloch equation with phase of pulse, without B1
        %and no slice correction. 
        disp( ['Corrected R2s estimation with Fn estimated from Bloch equations, no B1, no slice correction ', file_id,  '...']); 
        path=[path_results '/R2s_corr/BlochNoB1_NoSlcCorr/']; 
        if ~exist(path, 'dir')
                mkdir(path); 
        end

        R2s_corr_path = [path, 'R2s_BlochNoB1_noSlcCorr_', file_id, '.nii'];
        A_corr_path = [path, 'A_BlochNoB1_noSlcCorr_', file_id, '.nii'];
        res_norm_path = [path, 'resnorm_BlochNoB1_noSlcCorr_', file_id, '.nii'];
        residuals_path = [path, 'residuals_BlochNoB1_noSlcCorr_', file_id, '.nii'];
        
        %single slice option for faster calculation
        if opts.sgl_sag_slc == 1; 
        
           mask_opt = zeros(size(mag,1), size(mag,2), size(mag,3));
           mask_opt(:,size(mask,2)/2,:) = 1; 
           mask_opt = mask_opt.*double(mask); 
        else
           mask_opt = mask; 
        end

        if ~exist(R2s_corr_path, 'file')

            [R2s_corr_BlochNoB1NoSlcCorr, A_corr_BlochNoB1NoSlcCorr, resnorm_map_BlochNoB1NoSlcCorr, residuals_map_BlochNoB1NoSlcCorr] ...
                = CorrectedR2sMapEstimationFnNonLinFit(mag, Fn_BlochNoB1NoSlcCorr, te, mask_opt);

            mkdir(path); 
            nii = opts.nii_template; 

            nii.img = R2s_corr_BlochNoB1NoSlcCorr; 
            save_untouch_nii(nii, R2s_corr_path);
            nii.img =  A_corr_BlochNoB1NoSlcCorr; 
            save_untouch_nii(nii, A_corr_path);
            nii.img =  resnorm_map_BlochNoB1NoSlcCorr; 
            save_untouch_nii(nii, res_norm_path);
            nii.img =  residuals_map_BlochNoB1NoSlcCorr; 
            save_untouch_nii(nii, residuals_path);

        else
            tmp = load_untouch_nii(R2s_corr_path); 
            R2s_corr_BlochNoB1NoSlcCorr = double(tmp.img); 
            tmp = load_untouch_nii(A_corr_path); 
            A_corr_BlochNoB1NoSlcCorr = double(tmp.img); 
         %   tmp = load_untouch_nii(res_norm_path); 
          %  resnorm_map_BlochNoB1NoSlcCorr = double(tmp.img); 
        end

        res.R2s_corr_BlochNoB1NoSlcCorr = R2s_corr_BlochNoB1NoSlcCorr; 
        res.A_corr_BlochNoB1NoSlcCorr = A_corr_BlochNoB1NoSlcCorr; 
      %  res.resnorm_map_BlochNoB1NoSlcCorr = resnorm_map_BlochNoB1NoSlcCorr; 
      

%     
%     %% Estimate a corrected R2s with Fn from pulse shape wit Preibisch model
%     
%     disp( ['Corrected R2s estimation with Preibisch model ', file_id,  '...']); 
%     path=[path_results '/R2s_corr/Pulse_shape_Preibisch/']; 
%     if ~exist(path, 'dir')
%             mkdir(path); 
%     end
%     
%     R2s_corr_path = [path, 'R2s_corr_pulseshape_Preibisch_', file_id, '.nii'];
%     A_corr_path = [path, 'A_pulseshape_Preibisch_', file_id, '.nii'];
%     res_norm_path = [path, 'resnorm_pulseshape_Preibisch_', file_id, '.nii'];
%     residuals_path = [path, 'residuals_pulseshape_Preibisch_', file_id, '.nii'];
%     
%     if ~exist(R2s_corr_path, 'file')
%         
%         [R2s_corr_Preibisch, A_corr_Preibisch, resnorm_map_corr_Preibisch, residuals_map_corr_Preibisch] ...
%             = CorrectedR2sMapEstimationFnNonLinFit(mag, Fn_preibisch, te, mask);
% 
%         mkdir(path); 
%         nii = opts.nii_template; 
% 
%         nii.img = R2s_corr_Preibisch; 
%         save_untouch_nii(nii, R2s_corr_path);
%         nii.img =  A_corr_Preibisch; 
%         save_untouch_nii(nii, A_corr_path);
%         nii.img =  resnorm_map_corr_Preibisch; 
%         save_untouch_nii(nii, res_norm_path);
%         nii.img =  residuals_map_corr_Preibisch; 
%         save_untouch_nii(nii, residuals_path);
%         
%     else
%         tmp = load_untouch_nii(R2s_corr_path); 
%         R2s_corr_Preibisch = double(tmp.img); 
%         tmp = load_untouch_nii(A_corr_path); 
%         A_corr_Preibisch = double(tmp.img); 
%         tmp = load_untouch_nii(res_norm_path); 
%         resnorm_map_corr_Preibisch = double(tmp.img); 
%     end
%     
%     res.R2s_corr_Preibisch = R2s_corr_Preibisch; 
%     %res.A_corr_Preibisch = R2s_corr_Preibisch; 
%     %res.resnorm_map_corr_Preibisch = resnorm_map_corr_Preibisch; 
%     
end

