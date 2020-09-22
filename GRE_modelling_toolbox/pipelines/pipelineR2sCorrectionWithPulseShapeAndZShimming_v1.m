function [ res ] = pipelineR2sCorrectionWithPulseShapeAndZShimming_v1(mag, phase, te, MzShim, pulse, file_id, path_results, path_src_nii, z0, opts, B1_map)
%PIPELINER2SCORRECTIONWITHPULSESHAPE Pipeline for R2s correction for
%z-shimmed data. 
%   @param   MzShim          Table containing zShim Moments for each slice
%                            and echo [Nslcie, Nte] [µs mT/m]

    gamma = 267.51;                      % gyromagnetic ratio 42.57*2*pi in MHz/T
    
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

              opts_bet = '-m -f 0.35 -R -S -B -Z';
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
    % the phase is already unrwapped with PRELUDE!
    
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
                  ' -o ',  unwrapped, ' -s ', ' -m ', mask_path]; 
              system(prelude_call); 

              %load unwrapped mask 
              tmp = load_untouch_nii(unwrapped); 
              phase_prelude = double(tmp.img); 
          else 
              tmp = load_untouch_nii(unwrapped); 
              phase_prelude = double(tmp.img); 
          end
    end

    res.phase_prelude = phase_prelude; 

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
                te_selected = te(1:6); 

                bUnwrap = 1; 
                [ dw0, phi0] = LinearFitOfB0map(phase,mask, te_selected, 0, bUnwrap);

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
            if ~exist(path_B0_12, 'file')
          
                dphi =  phase_prelude(:,:,:,2) - phase_prelude(:,:,:,1);
                dTE = (te(2) - te(1))*1E-3; 
                dw0_12 = dphi./dTE; 

                nii = opts.nii_template; 
                nii.img =  dw0_12; 
                save_untouch_nii(nii, path_B0_12);
            else
                tmp = load_untouch_nii([path_B0_12]); 
                dw0_12 = double(tmp.img); 
            end
            
            path=[path_results, '/B0_map/B0_from_prelude/Echo13/'];
            if ~exist(path, 'dir')
                mkdir(path); 
            end
            
            path_B0_13 = [path, 'B0_from_prelude_13_', file_id, '.nii'];
            if ~exist(path_B0_13, 'file')
                dTE =  (te(3) -  te(1))*1E-3;
                dphi = phase_prelude(:,:,:,3) - phase_prelude(:,:,:,1);
                dw0_13 = dphi./dTE; 
                
                nii = opts.nii_template; 
                nii.img =  dw0_13; 
                save_untouch_nii(nii, path_B0_13);
            else
                tmp = load_untouch_nii(path_B0_13); 
                dw0_13 = double(tmp.img); 
            end
            
            dw0 = dw0_12; 
            res.dw0_12 = dw0_12; 
            res.dw0_13 = dw0_13;
   end
   

    %% Do a linear fit of the first four echo in order to estimate phi_0 
    %of the phase signal 

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
        [ dw0_lin, phi0_lin] = LinearFitOfB0map(phase_prelude,mask, te_selected, 0, bUnwrap);

        nii = opts.nii_template;
        nii.img = phi0_lin; 
        save_untouch_nii(nii, phi0_path); 
        nii.img = dw0_lin; 
        save_untouch_nii(nii, dw0_path); 
     else
       tmp = load_untouch_nii(phi0_path); 
       phi0_lin = double(tmp.img); 
       tmp = load_untouch_nii(dw0_path); 
       dw0_lin = double(tmp.img); 
    end
    res.dw0_lin = dw0_lin; 
    res.phi0_lin = phi0_lin; 

      
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

    %% 
    % ----------------------------------------------------------------
    % Estiamte Fn from z-shimmed data without slice correction and 
    % without B1 correction
    % ----------------------------------------------------------------

        
    disp( ['Estimation of Fn from z-shimmed data ', file_id,  '...']); 
    path=[path_results '/Fn_maps/zShimmed_Bloch_noB1_noSlcCorr/']; 
    if ~exist(path, 'dir')
            mkdir(path); 
    end

    Fn_map_bloch_path = [path, 'Fn_map_zShimmed_Bloch_noB1_noSlcCorr_', opts.B0_method,'_', file_id , '.nii'];

    if ~exist(Fn_map_bloch_path, 'file')
        if bB1_map == 1
            bSlcCorr = 0;
        else
            bSlcCorr = 1; %no slice correction due to the field gradient
        end
        Gsus = Gz/(gamma*z0); %convert from rad/(s*slice_thicness) to mT/m
        [ Fn_zShim_noBl_noSlcCorr] =  estimateFnMapFromCplxSliceProfileWithBlochEqnAndZShimming_v1(pulse, Gsus, MzShim, z0,  mask, te, bSlcCorr);

        nii = opts.nii_template_4D;  
        nii.img = Fn_zShim_noBl_noSlcCorr; 
        save_untouch_nii(nii, Fn_map_bloch_path);
    else
        tmp = load_untouch_nii(Fn_map_bloch_path); 
        Fn_zShim_noBl_noSlcCorr = double(tmp.img); 

    end
    res.Fn_zShim_noBl_noSlcCorr = Fn_zShim_noBl_noSlcCorr;
   
    
    %% 
    % ----------------------------------------------------------------
    % Estiamte corrected R2s with Fn from z-shimmed data without slice 
    % correction and without B1 correction
    % ----------------------------------------------------------------
    
    
    disp( ['Corrected R2s estimation with Fn estimated from zShimming, no B1 and no slice correction ', file_id,  '...']); 
    path=[path_results '/R2s_corr/BlochNoB1_noSlcCorr/']; 
    if ~exist(path, 'dir')
            mkdir(path); 
    end

    R2s_corr_path = [path, 'R2s_zShim_Bloch_noB1_noSlcCorr_', file_id, '.nii'];
    A_corr_path = [path, 'A_zShim_Bloch_noB1_noSlcCorr_', file_id, '.nii'];
    res_norm_path = [path, 'resnorm_zShim_Bloch_noB1_noSlcCorr_', file_id, '.nii'];
    residuals_path = [path, 'residuals_zShim_Bloch_noB1_noSlcCorr_', file_id, '.nii'];

    if ~exist(R2s_corr_path, 'file')

        [R2s_corr_zShimNoB1NoSlcCorr, A_corr_zShimNoB1NoSlcCorr, resnorm_map_zShimNoB1NoSlcCorr, residuals_map_zShimNoB1NoSlcCorr] ...
            = CorrectedR2sMapEstimationFnNonLinFit(mag, Fn_zShim_noBl_noSlcCorr, te, mask);

        mkdir(path); 
        nii = opts.nii_template; 

        nii.img = R2s_corr_zShimNoB1NoSlcCorr; 
        save_untouch_nii(nii, R2s_corr_path);
        nii.img =  A_corr_zShimNoB1NoSlcCorr; 
        save_untouch_nii(nii, A_corr_path);
        nii.img =  resnorm_map_zShimNoB1NoSlcCorr; 
        save_untouch_nii(nii, res_norm_path);
        nii.img =  residuals_map_zShimNoB1NoSlcCorr; 
        save_untouch_nii(nii, residuals_path);

    else
        tmp = load_untouch_nii(R2s_corr_path); 
        R2s_corr_zShimNoB1NoSlcCorr = double(tmp.img); 
        tmp = load_untouch_nii(A_corr_path); 
        A_corr_zShimNoB1NoSlcCorr = double(tmp.img); 
        tmp = load_untouch_nii(res_norm_path); 
        resnorm_map_zShimNoB1NoSlcCorr = double(tmp.img); 
    end

    res.R2s_corr_zShimNoB1NoSlcCorr = R2s_corr_zShimNoB1NoSlcCorr; 
    res.A_corr_zShimNoB1NoSlcCorr = A_corr_zShimNoB1NoSlcCorr; 
    res.resnorm_map_zShimNoB1NoSlcCorr = resnorm_map_zShimNoB1NoSlcCorr; 
    
    
    
    if bB1_map
        
        %% 
        % ----------------------------------------------------------------
        % Estiamte Fn from z-shimmed data without slice correction and 
        % with B1 correction
        % ----------------------------------------------------------------


        disp( ['Estimation of Fn from z-shimmed data with B1 without slc corr ', file_id,  '...']); 
        path=[path_results '/Fn_maps/zShimmed_Bloch_withB1_noSlcCorr/']; 
        if ~exist(path, 'dir')
                mkdir(path); 
        end

        Fn_map_bloch_path = [path, 'Fn_map_zShimmed_Bloch_withB1_noSlcCorr_', opts.B0_method,'_', file_id , '.nii'];

        if ~exist(Fn_map_bloch_path, 'file')
   
            bSlcCorr = 0;
       
            Gsus = Gz/(gamma*z0); %convert from rad/(s*slice_thicness) to mT/m
            [ Fn_zShim_withB1_noSlcCorr] =  estimateFnMapFromCplxSliceProfileWithBlochEqnAndZShimming_v1(pulse, Gsus, MzShim, z0,  mask, te, bSlcCorr, B1_map);

            nii = opts.nii_template_4D;  
            nii.img = Fn_zShim_withB1_noSlcCorr; 
            save_untouch_nii(nii, Fn_map_bloch_path);
        else
            tmp = load_untouch_nii(Fn_map_bloch_path); 
            Fn_zShim_withB1_noSlcCorr = double(tmp.img); 

        end
        res.Fn_zShim_withBl_noSlcCorr = Fn_zShim_withB1_noSlcCorr;    



        %% 
        % ----------------------------------------------------------------
        % Estiamte corrected R2s with Fn from z-shimmed data without slice 
        % correction and without B1 correction
        % ----------------------------------------------------------------


        disp( ['Corrected R2s estimation with Fn estimated from zShimming, no B1 and no slice correction ', file_id,  '...']); 
        path=[path_results '/R2s_corr/BlochWithB1_noSlcCorr/']; 
        if ~exist(path, 'dir')
                mkdir(path); 
        end

        R2s_corr_path = [path, 'R2s_zShim_Bloch_withB1_noSlcCorr_', file_id, '.nii'];
        A_corr_path = [path, 'A_zShim_Bloch_withB1_noSlcCorr_', file_id, '.nii'];
        res_norm_path = [path, 'resnorm_zShim_Bloch_withB1_noSlcCorr_', file_id, '.nii'];
        residuals_path = [path, 'residuals_zShim_Bloch_withB1_noSlcCorr_', file_id, '.nii'];

        if ~exist(R2s_corr_path, 'file')

            [R2s_corr_zShimWithB1NoSlcCorr, A_corr_zShimWithB1NoSlcCorr, resnorm_map_zShimwithB1NoSlcCorr, residuals_map_zShimWithB1NoSlcCorr] ...
                = CorrectedR2sMapEstimationFnNonLinFit(mag, Fn_zShim_withB1_noSlcCorr, te, mask);

            mkdir(path); 
            nii = opts.nii_template; 

            nii.img = R2s_corr_zShimWithB1NoSlcCorr; 
            save_untouch_nii(nii, R2s_corr_path);
            nii.img =  A_corr_zShimWithB1NoSlcCorr; 
            save_untouch_nii(nii, A_corr_path);
            nii.img =  resnorm_map_zShimwithB1NoSlcCorr; 
            save_untouch_nii(nii, res_norm_path);
            nii.img =  residuals_map_zShimWithB1NoSlcCorr; 
            save_untouch_nii(nii, residuals_path);

        else
            tmp = load_untouch_nii(R2s_corr_path); 
            R2s_corr_zShimWithB1NoSlcCorr = double(tmp.img); 
            tmp = load_untouch_nii(A_corr_path); 
            A_corr_zShimWithB1NoSlcCorr = double(tmp.img); 
            tmp = load_untouch_nii(res_norm_path); 
            resnorm_map_zShimwithB1NoSlcCorr = double(tmp.img); 
        end

        res.R2s_corr_zShimWithB1NoSlcCorr = R2s_corr_zShimWithB1NoSlcCorr; 
        res.A_corr_zShimNoB1NoSlcCorr = A_corr_zShimWithB1NoSlcCorr; 
        res.resnorm_map_zShimNoB1WithSlcCorr = resnorm_map_zShimwithB1NoSlcCorr; 


        
        
        
        
        
        
        
        %% 
        % ----------------------------------------------------------------
        % Estiamte Fn from z-shimmed data with slice correction and 
        % with B1 correction
        % ----------------------------------------------------------------


        disp( ['Estimation of Fn from z-shimmed data with B1 and slc corr ', file_id,  '...']); 
        path=[path_results '/Fn_maps/zShimmed_Bloch_withB1_withSlcCorr/']; 
        if ~exist(path, 'dir')
                mkdir(path); 
        end

        Fn_map_bloch_path = [path, 'Fn_map_zShimmed_Bloch_withB1_withSlcCorr_', opts.B0_method,'_', file_id , '.nii'];

        if ~exist(Fn_map_bloch_path, 'file')
   
            bSlcCorr = 1;
       
            Gsus = Gz/(gamma*z0); %convert from rad/(s*slice_thicness) to mT/m
            [ Fn_zShim_withB1_withSlcCorr] =  estimateFnMapFromCplxSliceProfileWithBlochEqnAndZShimming_v1(pulse, Gsus, MzShim, z0,  mask, te, bSlcCorr, B1_map);

            nii = opts.nii_template_4D;  
            nii.img = Fn_zShim_withB1_withSlcCorr; 
            save_untouch_nii(nii, Fn_map_bloch_path);
        else
            tmp = load_untouch_nii(Fn_map_bloch_path); 
            Fn_zShim_withB1_withSlcCorr = double(tmp.img); 

        end
        res.Fn_zShim_withB1_withSlcCorr = Fn_zShim_withB1_withSlcCorr;    



        %% 
        % ----------------------------------------------------------------
        % Estiamte corrected R2s with Fn from z-shimmed data without slice 
        % correction and without B1 correction
        % ----------------------------------------------------------------


        disp( ['Corrected R2s estimation with Fn estimated from zShimming, with B1 and with slice correction ', file_id,  '...']); 
        path=[path_results '/R2s_corr/BlochWithB1_withSlcCorr/']; 
        if ~exist(path, 'dir')
                mkdir(path); 
        end

        R2s_corr_path = [path, 'R2s_zShim_Bloch_withB1_withSlcCorr_', file_id, '.nii'];
        A_corr_path = [path, 'A_zShim_Bloch_withB1_withSlcCorr_', file_id, '.nii'];
        res_norm_path = [path, 'resnorm_zShim_Bloch_withB1_withSlcCorr_', file_id, '.nii'];
        residuals_path = [path, 'residuals_zShim_Bloch_withB1_withSlcCorr_', file_id, '.nii'];

        if ~exist(R2s_corr_path, 'file')

            [R2s_corr_zShimWithB1WithSlcCorr, A_corr_zShimWithB1WithSlcCorr, resnorm_map_zShimWithB1WithSlcCorr, residuals_map_zShimWithB1WithSlcCorr] ...
                = CorrectedR2sMapEstimationFnNonLinFit(mag, Fn_zShim_withB1_withSlcCorr, te, mask);

            mkdir(path); 
            nii = opts.nii_template; 

            nii.img = R2s_corr_zShimWithB1WithSlcCorr; 
            save_untouch_nii(nii, R2s_corr_path);
            nii.img =  A_corr_zShimWithB1WithSlcCorr; 
            save_untouch_nii(nii, A_corr_path);
            nii.img =  resnorm_map_zShimWithB1WithSlcCorr; 
            save_untouch_nii(nii, res_norm_path);
            nii.img =  residuals_map_zShimWithB1WithSlcCorr; 
            save_untouch_nii(nii, residuals_path);

        else
            tmp = load_untouch_nii(R2s_corr_path); 
            R2s_corr_zShimWithB1WithSlcCorr = double(tmp.img); 
            tmp = load_untouch_nii(A_corr_path); 
            A_corr_zShimWithB1WithSlcCorr = double(tmp.img); 
            tmp = load_untouch_nii(res_norm_path); 
            resnorm_map_zShimWithB1WithSlcCorr = double(tmp.img); 
        end

        res.R2s_corr_zShimWithB1WithSlcCorr = R2s_corr_zShimWithB1WithSlcCorr; 
        res.A_corr_zShimNoB1NoSlcCorr = A_corr_zShimWithB1WithSlcCorr; 
        res.resnorm_map_zShimWithB1WithSlcCorr = resnorm_map_zShimWithB1WithSlcCorr; 


        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    end
    
    
    
%         
% %     
% %    %% 
% %     % ----------------------------------------------------------------
% %     % Estiamte corrected R2s with Fn from z-shimmed from the choes where no
% %     % z-shimming was applied 
% %     % ----------------------------------------------------------------
% %     
% %     
% %         
% %     disp( ['Corrected R2s estimation with Fn estimated only from echoes estimaed without zShim, no B1 and no slice correction ', file_id,  '...']); 
% %     path=[path_results '/R2s_corr/NoZShim_FirstEchoes_BlochNoB1_wtihSlcCorr/']; 
% %     if ~exist(path, 'dir')
% %             mkdir(path); 
% %     end
% % 
% %     R2s_corr_path = [path, 'R2s_zShim_firstEchoes_Bloch_noB1_noSlcCorr_', file_id, '.nii'];
% %     A_corr_path = [path, 'A_zShim_firstEchoes_Bloch_noB1_noSlcCorr_', file_id, '.nii'];
% %     res_norm_path = [path, 'resnorm_zShim_firstEchoes_Bloch_noB1_noSlcCorr_', file_id, '.nii'];
% %     residuals_path = [path, 'residuals_zShim_firstEchoes_Bloch_noB1_noSlcCorr_', file_id, '.nii'];
% % 
% %     if ~exist(R2s_corr_path, 'file')
% % 
% %         startEcho = opts.opts_zShim.startEcho;
% %         if (opts.opts_zShim.pattern_nr == -1) 
% %             te_fit = te;          
% %             mag_fit = mag; 
% %             Fn_zShim_noBl_noSlcCorr_fit = Fn_zShim_noBl_noSlcCorr; 
% %             mask_fit = mask; 
% %             
% %         else
% %             te_fig = te(1:startEcho);
% %             mag_fit = mag(:,:,1:startEcho); 
% %             Fn_zShim_noBl_noSlcCorr_fit = Fn_zShim_noBl_noSlcCorr(:,:,1:startEcho); 
% %             mask_fit = mask(:,:,1:startEcho); 
% %         end
% %         
% %         [R2s_corr_firstEchoes_zShimNoB1NoSlcCorr, A_corr_firstEchoes_zShimNoB1NoSlcCorr, resnorm_map_firstEchoes_zShimNoB1NoSlcCorr, residuals_map_firstEchoes_zShimNoB1NoSlcCorr] ...
% %             = CorrectedR2sMapEstimationFnNonLinFit(mag_fit, Fn_zShim_noBl_noSlcCorr_fit, te_fit, mask_fit);
% % 
% %         mkdir(path); 
% %         nii = opts.nii_template; 
% % 
% %         nii.img = R2s_corr_firstEchoes_zShimNoB1NoSlcCorr; 
% %         save_untouch_nii(nii, R2s_corr_path);
% %         nii.img =  A_corr_firstEchoes_zShimNoB1NoSlcCorr; 
% %         save_untouch_nii(nii, A_corr_path);
% %         nii.img =  resnorm_map_firstEchoes_zShimNoB1NoSlcCorr; 
% %         save_untouch_nii(nii, res_norm_path);
% %         nii.img =  residuals_map_firstEchoes_zShimNoB1NoSlcCorr; 
% %         save_untouch_nii(nii, residuals_path);
% % 
% %     else
% %         tmp = load_untouch_nii(R2s_corr_path); 
% %         R2s_corr_firstEchoes_zShimNoB1NoSlcCorr = double(tmp.img); 
% %         tmp = load_untouch_nii(A_corr_path); 
% %         A_corr_firstEchoes_zShimNoB1NoSlcCorr = double(tmp.img); 
% %         tmp = load_untouch_nii(res_norm_path); 
% %         resnorm_map_firstEchoes_zShimNoB1NoSlcCorr = double(tmp.img); 
% %     end
% % 
% %     res.R2s_corr_firstEchoes_zShimNoB1NoSlcCorr = R2s_corr_firstEchoes_zShimNoB1NoSlcCorr; 
% %     res.A_corr_firstEchoes_zShimNoB1NoSlcCorr = A_corr_firstEchoes_zShimNoB1NoSlcCorr; 
% %     res.resnorm_map_firstEchoes_zShimNoB1NoSlcCorr = resnorm_map_firstEchoes_zShimNoB1NoSlcCorr; 
% %     
% %     
% %     
% %  

% 
%     %% 
%     % ----------------------------------------------------------------
%     % Estiamte R2s with CUDA fit only in case that no z-shim is applied 
%     % (refrence acquistion without z-shim gradients)
%     % ----------------------------------------------------------------
%     
%     addpath('/media/data/physics/lukas/MATLAB/MRIDOG/CANDI/tin_R2star/matlab/clmf_expr2/');
%     disp(['Estimate R2s with CUDA ', file_id,  '...']); 
% 
%     GRE = path_src_nii; 
%     path=[path_results, '/R2s_cuda_cmd/']; 
%     if ~exist(path); 
%         mkdir(path); 
%     end
% 
%     startEcho = opts.opts_zShim.startEcho;
% 
%     te_cuda = te(1:startEcho);
% 
%     
%     
%     dst_nii = [path, 'R2s_cuda_', num2str(length(te_cuda)), 'E_', file_id, '.nii.gz'];
%     if ~exist(dst_nii, 'file'); 
%               
%         
%         %call "pure" Cuda fit without any preprocessing/filtering
%         callCUDAR2sFit(GRE, dst_nii, te_cuda);
%         
%         %load map
%         tmp = load_untouch_nii(dst_nii); 
%         R2s_cuda = double(tmp.img); 
%         R2s_cuda(isinf(R2s_cuda)) = 0;
%         R2s_cuda(isnan(R2s_cuda)) = 0;
%         
%     else 
%         tmp = load_untouch_nii(dst_nii); 
%         R2s_cuda= double(tmp.img); 
%         R2s_cuda(isinf(R2s_cuda)) = 0;
%         R2s_cuda(isnan(R2s_cuda)) = 0;
%     end
%     
%     res.R2s_cuda = R2s_cuda; 

%     
%     
% 
%     
%     

end

