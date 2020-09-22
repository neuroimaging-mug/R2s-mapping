function [ res ] = estimateGradientMaps( mag, phase, te, file_id, path_results, opts )
%ESTIMATEGRADIENTMAPS Summary of this function goes here
%   Detailed explanation goes here


    gamma = 267.51;                      % gyromagnetic ratio 42.57*2*pi in MHz/T
    

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
    used_mask_path = []; 
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
         used_mask_path = bet_mask_path; 
         tmp = load_untouch_nii(bet_mask_path); 
         mask = double(tmp.img); 
         
    else % Create simple mask 
            path=[path_results, '/Mask/']; 
            if exist(path, 'dir') == 0
                mkdir(path);   
            end
            mask_path = [path, 'Mask_', file_id, '.nii'];
            used_mask_path = mask_path; 
            if exist(mask_path, 'file') == 0
                %create a simple mask 
                mask_tmp =  mag(:,:,:,1);  
                mask_tmp(mask_tmp < 100) = 0;  
                mask_tmp(mask_tmp >= 100) = 1; 
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
                  ' -o ',  unwrapped, ' -s', ' -m ', used_mask_path]; 
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
        
        mask_mag =  mag; 
        mask_mag(mask_mag <= 75) = 0; 
        mask_mag(mask_mag > 75) = 1; 
        mask_mag = mask_mag(:,:,:,1).*mask; 

        
        SE = strel('disk',4);
        mask_er = imerode(mask_mag, SE); 
        [ Gx, Gy, ~] = gradient(dw0);
        Gx = Gx.*mask_er; 
        Gy = Gy.*mask_er; 
        
        

        Gz =  zeros(size(dw0)); 
        nSlices = size(dw0,3); 


        for i=1:nSlices
            Gtmp = zeros(size(dw0,1), size(dw0,2)); 
            if i == 1
               Gz(:,:,1) =  dw0(:,:,2) - dw0(:,:,1); 
            elseif i== nSlices
               Gz(:,:,nSlices) =  dw0(:,:,nSlices) - dw0(:,:,nSlices-1); 
            else

               d_central = 0.5*(dw0(:,:,i+1) - dw0(:,:,i-1)); 
               d_pos = dw0(:,:,i+1) - dw0(:,:,i); 
               d_neg = dw0(:,:,i) - dw0(:,:,i-1); 

               mask_sum = sum(mask_er(:,:,i-1:i+1),3); 

               mask_neg =  mask_er(:,:,i-1) - mask_er(:,:,i)  -  mask_er(:,:,i+1); 
               mask_pos =  mask_er(:,:,i+1) - mask_er(:,:,i)  -  mask_er(:,:,i-1); 

               Gtmp(mask_sum == 3) = d_central(mask_sum == 3); 
%                Gtmp(mask_pos == 0) = d_pos(mask_pos == 0); 
%                Gtmp(mask_neg == 0) = d_neg(mask_neg == 0); 


               Gz(:,:,i) = Gtmp; 
            end

        end

        
        
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

end

