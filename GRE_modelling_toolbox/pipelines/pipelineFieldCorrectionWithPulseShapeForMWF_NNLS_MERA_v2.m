function [ res ] = pipelineFieldCorrectionWithPulseShapeForMWF_NNLS_MERA_v2(mag, phase, te, pulse, file_id, path_results, path_src_nii, z0, opts, B1_map)
%PIPELINER2SCORRECTIONWITHPULSESHAPE Pipeline for correction of field
%influences for MWF-mapping. 
%v2: Regualrizatin paramter beta_reg for NNLS is added as addtional
%paramter. 
    gamma = 267.51;                      % gyromagnetic ratio 42.57*2*pi in MHz/T
    %Check if B1-map is available 
    bB1_map = 0; 
    if nargin == 10
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
         bet_mask_path = [path, 'bet_mag_echo1_', file_id, '_mask.nii.gz']; 
         tmp = load_untouch_nii(bet_mask_path); 
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
                mask_tmp(mask_tmp < 50) = 0;  
                mask_tmp(mask_tmp >= 50) = 1; 
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

              %mask_path = [path_results, '/Mask/', 'mask_', file_id]; 

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
            path_phi0_12 = [path, 'phi0_from_prelude_12_', file_id, '.nii'];
            if ~exist(path_B0_12, 'file') || ~exist(path_phi0_12, 'file')
          
                dTE = (te(2) - te(1))*1E-3; 
                if opts.bprelude == 1
                    dphi = phase_prelude(:,:,:,2) - phase_prelude(:,:,:,1);
                    dw0_12 = dphi./dTE; 
                    phi0_12 = phase_prelude(:,:,:,2) - dw0_12*te(2)*1E-3;
                else
                    dphi = phase(:,:,:,2) - phase(:,:,:,1);
                    dw0_12 = dphi./dTE; 
                    phi0_12 = phase(:,:,:,2) - dw0_12*te(2)*1E-3;
                end
       
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
                
                if opts.bprelude == 1
                    dphi = phase_prelude(:,:,:,3) - phase_prelude(:,:,:,1);
                    dw0_13 = dphi./dTE; 
                    phi0_13 = phase_prelude(:,:,:,2) - dw0_13*te(2)*1E-3;
                else
                    dphi = phase(:,:,:,3) - phase(:,:,:,1);
                    dw0_13 = dphi./dTE; 
                    phi0_13 = phase(:,:,:,2) - dw0_13*te(2)*1E-3;
                end
       
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
            
            dw0 = dw0_13; 
            phi0 = phi0_13; 
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
    % Estimate MWF with MERA NNLS
    % ----------------------------------------------------------------

%     mask_tmp = zeros(size(mask)); 
%     mask_tmp(128:138, 128, :) = 1; 

  
    disp( ['MWF estimation from NNLS fit ', file_id,  '...']); 
    
    reg_str = num2str(opts.beta_reg); reg_str(reg_str == '.') = 'd'; 
    path=[path_results '/MWF/NNLS_MERA_reg_', reg_str, '/MWF_MERA/']; 
    if ~exist(path, 'dir')
            mkdir(path); 
    end

    M0_path = [path, 'M0_MERA_', file_id, '.nii'];
    MWF_path = [path, 'MWF_MERA_', file_id, '.nii'];

    if ~exist(MWF_path, 'file')

        [MWF_MERA_map,  M0]...
             = NNLS_MultiExponentialMapEstimation_MERA(mag, te, mask, opts.beta_reg);

        mkdir(path); 
        nii = opts.nii_template; 
        
        nii.img = M0; 
        save_untouch_nii(nii, M0_path); 

        nii = opts.nii_template; 
        nii.img = MWF_MERA_map; 
        save_untouch_nii(nii, MWF_path);

    else
        tmp = load_untouch_nii(M0_path); 
        M0 = double(tmp.img); 
        tmp = load_untouch_nii(MWF_path); 
        MWF_MERA_map = double(tmp.img);             
    end

    res.M0 = M0; 
    res.MWF_MERA_map = MWF_MERA_map;


    
    
    
    
    

    %     ----------------------------------------------------------------
    %     Estimate Fn with Bloch equations and only the mag of the slice 
    %     profile
    %     ----------------------------------------------------------------


    disp( ['Estimation of Fn by solving the Bloch equations for the magnitude slice profile ', file_id,  '...']); 
    path=[path_results '/Fn_maps/BlochMagOnly/']; 
    if ~exist(path, 'dir')
            mkdir(path); 
    end

    Fn_map_bloch_path = [path, 'Fn_map_BlochMagOnly_', opts.B0_method,'_', file_id , '.nii'];

    if ~exist(Fn_map_bloch_path, 'file')

        Gsus = Gz/(gamma*z0); %convert from rad/(s*slice_thicness) to mT/m
        [ Fn_BlochMagOnly] = estimateFnMapFromMagSliceProfileWithBlochEqn(pulse, Gsus, z0,  mask, te);

        nii = opts.nii_template_4D;  
        nii.img = Fn_BlochMagOnly; 
        save_untouch_nii(nii, Fn_map_bloch_path);
    else
        tmp = load_untouch_nii(Fn_map_bloch_path); 
        Fn_BlochMagOnly = double(tmp.img); 

    end
    res.Fn_BlochMagOnly = Fn_BlochMagOnly;
    
    
    
    
%     
%     % ----------------------------------------------------------------
%     % Estimate a corrected MWF with MERA NNLS using magnitude only model
%     % ----------------------------------------------------------------
%     
%     
% %     mask_tmp = zeros(size(mask)); 
% %     mask_tmp(100:120,100:120,10) = 1; 
%     
%     disp( ['MWF estimation from corrected data from Fmag without B1 and slc corr with  NNLS fit ', file_id,  '...']); 
%     
%     e_inv =[0.05, 0.08, 0.10:0.05:0.60]; 
%     e_inv = 0.60; 
%     for i=1:length(e_inv);  
%         
%         e_str = num2str(e_inv(i)); e_str(e_str == '.') = '_'; 
%         
% 
%         path=[path_results '/MWF/NNLS_MERA_reg_', reg_str, '/MWF_MERA_corr_Fmag/']; 
%         
%         if ~exist(path, 'dir')
%                 mkdir(path); 
%         end
% 
%         M0_path = [path, 'M0_MERA_corr_Fmag_', e_str, '_', file_id, '.nii'];
%         MWF_path = [path, 'MWF_MERA_corr_Fmag_', e_str, '_', file_id, '.nii'];
%         mask_corr_path = [path, 'mask_signal_corrr_Fmag_', e_str, '_', file_id, '.nii'];
%         if ~exist(MWF_path, 'file')
% 
%             [MWF_MERA_map_corr_Fmag{i}, M0_corr_Fmag{i}, mask_signal_Fmag{i}] = ...
%                 NNLS_MultiExponentialMapEstimation_MERA_withCorrection(mag, te, mask, Fn_BlochMagOnly, e_inv(i), opts.beta_reg);
% 
%             mkdir(path); 
%             nii = opts.nii_template; 
% 
%             nii.img = M0_corr_Fmag{i}; 
%             save_untouch_nii(nii, M0_path); 
% 
%             nii = opts.nii_template; 
%             nii.img = MWF_MERA_map_corr_Fmag{i}; 
%             save_untouch_nii(nii, MWF_path);
% 
%             mkdir(path); 
%             nii = opts.nii_template_4D; 
%             nii.img = mask_signal_Fmag{i}; 
%             save_untouch_nii(nii, mask_corr_path);
%         else
%             tmp = load_untouch_nii(M0_path); 
%             M0_corr_Fmag{i} = double(tmp.img); 
%             tmp = load_untouch_nii(MWF_path); 
%             MWF_MERA_map_corr_Fmag{i} = double(tmp.img);             
%             tmp = load_untouch_nii(mask_corr_path); 
%             mask_signal_Fmag{i} = double(tmp.img);           
%         end
% 
%     
%     end
% 
% 
%     res.M0_noB1_noSlcCorr = M0_corr_Fmag; 
%     res.MWF_MERA_noB1_noSlcCorr = MWF_MERA_map_corr_Fmag;
%     res.mask_noB1_noSlcCorr = mask_signal_Fmag; 
%  
%     
    
    
    
    
    

    % ----------------------------------------------------------------
    % Without B1 and without slice correction
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
        [ Fn_BlochNoB1NoSlcCorr] =  estimateFnMapFromCplxSliceProfileWithBlochEquations_v1(pulse, Gsus, z0,  mask, te, bSlcCorr);

        nii = opts.nii_template_4D;  
        nii.img = Fn_BlochNoB1NoSlcCorr; 
        save_untouch_nii(nii, Fn_map_bloch_path);
    else
        tmp = load_untouch_nii(Fn_map_bloch_path); 
        Fn_BlochNoB1NoSlcCorr = double(tmp.img); 

    end
    res.Fn_BlochNoB1NoSlcCorr = Fn_BlochNoB1NoSlcCorr;
    
    
   
    
    
    % ----------------------------------------------------------------
    % Estimate a corrected MWF with MERA NNLS
    % ----------------------------------------------------------------
    
    
%     mask_tmp = zeros(size(mask)); 
%     mask_tmp(100:120,100:120,10) = 1; 
    
    disp( ['MWF estimation from corrected data from Fn without B1 and slc corr with  NNLS fit ', file_id,  '...']); 
    
    e_inv =[0.05, 0.08, 0.10:0.05:0.60]; 
    e_inv = 0.60; 
    for i=1:length(e_inv);  
        
        e_str = num2str(e_inv(i)); e_str(e_str == '.') = '_'; 
        
        path=[path_results '/MWF/NNLS_MERA_reg_', reg_str, '/MWF_MERA_corr_Fn_noB1_noSlcCorr/']; 

        if ~exist(path, 'dir')
                mkdir(path); 
        end

        M0_path = [path, 'M0_MERA_corr_noB1_noSlcCorr_', e_str, '_', file_id, '.nii'];
        MWF_path = [path, 'MWF_MERA_corr_noB1_noSlcCorr_', e_str, '_', file_id, '.nii'];
        mask_corr_path = [path, 'mask_signal_corr_noB1_noSlcCorr_', e_str, '_', file_id, '.nii'];
        if ~exist(MWF_path, 'file')

            [MWF_MERA_map_corr{i}, M0_corr{i}, mask_signal{i}] = ...
                NNLS_MultiExponentialMapEstimation_MERA_withCorrection(mag, te, mask, Fn_BlochNoB1NoSlcCorr, e_inv(i), opts.beta_reg);

            mkdir(path); 
            nii = opts.nii_template; 

            nii.img = M0_corr{i}; 
            save_untouch_nii(nii, M0_path); 

            nii = opts.nii_template; 
            nii.img = MWF_MERA_map_corr{i}; 
            save_untouch_nii(nii, MWF_path);

            mkdir(path); 
            nii = opts.nii_template_4D; 
            nii.img = mask_signal{i}; 
            save_untouch_nii(nii, mask_corr_path);
        else
            tmp = load_untouch_nii(M0_path); 
            M0_corr{i} = double(tmp.img); 
            tmp = load_untouch_nii(MWF_path); 
            MWF_MERA_map_corr{i} = double(tmp.img);             
            tmp = load_untouch_nii(mask_corr_path); 
            mask_signal{i} = double(tmp.img);           
        end

    
    end


    res.M0_noB1_noSlcCorr = M0_corr; 
    res.MWF_MERA_noB1_noSlcCorr = MWF_MERA_map_corr;
    res.mask_noB1_noSlcCorr = mask_signal; 
 
    
        
    % ----------------------------------------------------------------
    % With B1 and with slice correction
    % ----------------------------------------------------------------

    
    
    if bB1_map == 1
        
        disp( ['Estimation of Fn by solving the Bloch equations for the slice profile including B1 errors and slice correction', file_id,  '...']); 
        path=[path_results '/Fn_maps/BlochWithB1_withSlcCorr/']; 
        if ~exist(path, 'dir')
                mkdir(path); 
        end

        Fn_map_bloch_path = [path, 'Fn_map_BlochWithB1_withSlcCorr_', opts.B0_method,'_', file_id , '.nii'];
        if ~exist(Fn_map_bloch_path, 'file')

            bSlcCorr = 1; %slice correction due to the field gradient
            Gsus = Gz/(gamma*z0); %convert from rad/(s*slice_thicness) to mT/m
            z0_vary = z0; 
            [ Fn_BlochWithB1WithSlcCorr] =  estimateFnMapFromCplxSliceProfileWithBlochEquations_v1(pulse, Gsus, z0_vary,  mask, te, bSlcCorr, B1_map);
            

            nii = opts.nii_template_4D;  
            nii.img = Fn_BlochWithB1WithSlcCorr; 
            save_untouch_nii(nii, Fn_map_bloch_path);
        else
            tmp = load_untouch_nii(Fn_map_bloch_path); 
            Fn_BlochWithB1WithSlcCorr = double(tmp.img); 

        end
        res.Fn_BlochWithB1WithSlcCorr = Fn_BlochWithB1WithSlcCorr;
        
        
        disp( ['MWF estimation from corrected data from Fn with B1 and slc corr with  NNLS fit ', file_id,  '...']); 

        for i=1:length(e_inv);  

            e_str = num2str(e_inv(i)); e_str(e_str == '.') = '_'; 

            path=[path_results '/MWF/NNLS_MERA_reg_', reg_str, '/MWF_MERA_corr_Fn_withB1_withSlcCorr/']; 
            if ~exist(path, 'dir')
                    mkdir(path); 
            end

            M0_path = [path, 'M0_MERA_corr_withB1_withSlcCorr_', e_str, '_', file_id, '.nii'];
            MWF_path = [path, 'MWF_MERA_corr_withB1_withSlcCorr_', e_str, '_', file_id, '.nii'];
            mask_corr_path = [path, 'mask_signal_withB1_withSlcCorr_', e_str, '_', file_id, '.nii'];
            if ~exist(MWF_path, 'file')

%                 %debugging mask! 
%                 disp('mask for debuging is activated!'); 
%                 mask(:,:,1:7) = 0; 
%                 mask(:,:,9:end) = 0; 
                
                [MWF_MERA_map_corr{i}, M0_corr{i}, mask_signal{i}] = ...
                    NNLS_MultiExponentialMapEstimation_MERA_withCorrection(mag, te, mask, Fn_BlochWithB1WithSlcCorr, e_inv(i), opts.beta_reg);

                mkdir(path); 
                nii = opts.nii_template; 

                nii.img = M0_corr{i}; 
                save_untouch_nii(nii, M0_path); 

                nii = opts.nii_template; 
                nii.img = MWF_MERA_map_corr{i}; 
                save_untouch_nii(nii, MWF_path);

                mkdir(path); 
                nii = opts.nii_template_4D; 
                nii.img = mask_signal{i}; 
                save_untouch_nii(nii, mask_corr_path);
            else
                tmp = load_untouch_nii(M0_path); 
                M0_corr{i} = double(tmp.img); 
                tmp = load_untouch_nii(MWF_path); 
                MWF_MERA_map_corr{i} = double(tmp.img);             
                tmp = load_untouch_nii(mask_corr_path); 
                mask_signal{i} = double(tmp.img);           
            end


        end


        res.M0_wB1_wSlcCorr = M0_corr; 
        res.MWF_MERA_wB1_wSlcCorr = MWF_MERA_map_corr;
        res.mask_wB1_wSlcCorr = mask_signal; 


    
    end
    

    
    % ----------------------------------------------------------------
    % Estimate Gz_meas with the method proposed by Alonso-Ortiz 
    % ----------------------------------------------------------------
    

    if opts.AlonsoCorr == 1


        disp( ['Estimate Gz fit with the Alonso-Ortiz approach ', file_id,  '...']); 
        path=[path_results '/Gz_fit/']; 
        if ~exist(path, 'dir')
                mkdir(path); 
        end

        dwz_fit_path = [path, 'dwz_fit_', opts.B0_method,'_', file_id , '.nii'];
        T2s_fit_path = [path, 'T2s_fit_', opts.B0_method,'_', file_id , '.nii'];
        S0_fit_path = [path, 'S0_fit_', opts.B0_method,'_', file_id , '.nii'];
        Fn_sinc_fit_path = [path, 'Fn_sinc_fit_', opts.B0_method,'_', file_id , '.nii'];
        if ~exist(dwz_fit_path, 'file')

    %         Gsus = Gz/(z0_cm); %convert to rad/(s*cm)
    %         voxel_dim = [1,1,z0]; %mm
           dwz_ms = Gz*1E-3; %convert to rad/ms
           %te index for IEW water
           te_idx = find(te > 5); 

            [T2s_fit, dw_fit] = CorrectedT2sMapEstimationSincNonLinFit(squeeze(mag(:,:,:,te_idx:end)), dwz_ms, te(te_idx:end), mask);

            [Ny, Nx, Nz] = size(dw_fit); 
            te_mat = repmat(permute(te, [1,4,3,2]),[Ny, Nx, Nz,1]); 
            dw_fit_mat = repmat(dw_fit, [1,1,1, length(te)]); 

            Fn_sinc_fit = sinc(dw_fit_mat.*te_mat/(2*pi));

            nii = opts.nii_template; 
            nii.img = dw_fit; 
            save_untouch_nii(nii, dwz_fit_path);
            nii.img = T2s_fit; 
            save_untouch_nii(nii, T2s_fit_path);
    %         nii.img = S0_fit; 
    %         save_untouch_nii(nii, S0_fit_path);
            nii = opts.nii_template_4D; 
            nii.img = Fn_sinc_fit; 
            save_untouch_nii(nii, Fn_sinc_fit_path);
        else
            tmp = load_untouch_nii(dwz_fit_path); 
            dw_fit = double(tmp.img); 
            tmp = load_untouch_nii(T2s_fit_path); 
            T2s_fit = double(tmp.img); 
    %         tmp = load_untouch_nii(S0_fit_path); 
    %         S0_fit = double(tmp.img); 

            tmp = load_untouch_nii(Fn_sinc_fit_path); 
            Fn_sinc_fit = double(tmp.img); 

        end
        res.dw_fit = dw_fit;
        res.T2s_fit = T2s_fit;
    %     res.S0_fit = S0_fit;
        res.Fn_sinc_fit = Fn_sinc_fit;


       disp( ['MWF estimation from corrected Fn data with  NNLS fit wit Alonso-Ortiz correction ', file_id,  '...']); 

       e_inv = 0.25; 
       for i=1:length(e_inv);  

            e_str = num2str(e_inv(i)); e_str(e_str == '.') = '_'; 

            path=[path_results '/MWF/NNLS_MERA_corr_AlonsoOrtiz_', e_str, '/']; 
            if ~exist(path, 'dir')
                    mkdir(path); 
            end

            M0_path = [path, 'M0_MERA_corr_AlonsoOrtiz_', e_str, '_', file_id, '.nii'];
            MWF_path = [path, 'MWF_MERA_corr_AlonsoOrtiz_', e_str, '_', file_id, '.nii'];
            mask_corr_path = [path, 'mask_signal_corr_AlonsoOrtiz_', e_str, '_', file_id, '.nii'];
            if ~exist(MWF_path, 'file')

                [MWF_MERA_map_corr{i}, M0_corr{i}, mask_signal{i}] = ...
                    NNLS_MultiExponentialMapEstimation_MERA_withCorrection(mag, te, mask, abs(Fn_sinc_fit), e_inv(i));

                mkdir(path); 
                nii = opts.nii_template; 

                nii.img = M0_corr{i}; 
                save_untouch_nii(nii, M0_path); 

                nii = opts.nii_template; 
                nii.img = MWF_MERA_map_corr{i}; 
                save_untouch_nii(nii, MWF_path);

                mkdir(path); 
                nii = opts.nii_template_4D; 
                nii.img = mask_signal{i}; 
                save_untouch_nii(nii, mask_corr_path);
            else
                tmp = load_untouch_nii(M0_path); 
                M0_corr{i} = double(tmp.img); 
                tmp = load_untouch_nii(MWF_path); 
                MWF_MERA_map_corr{i} = double(tmp.img);             
                tmp = load_untouch_nii(mask_corr_path); 
                mask_signal{i} = double(tmp.img);           
            end


        end


        res.M0_corr_AlonsoOrtiz = M0_corr; 
        res.MWF_MERA_map_corr_AlonsoOrtiz  = MWF_MERA_map_corr;
        res.mask_signal_AlonsoOrtiz  = mask_signal; 
    
    
    end


    

    
end

