function [res ...
    ] = averageNavigatorData( dat_path_c, nii_headers, dst_raw, dst_navi, bdetrend,bloadExample)
%AVERAGENAVIGATORDATA Scripts loads Siemens raw data an performs a
%correction of data with the navigator echo (last echo). 
%Afterwards images are registered with FSL Flirst and averaged. 
%   @param  dat_path_c       cell with path to the Siemens data file 
%   @param  nii_headers      cell containing nifti headers 
%   @param  dst_raw          path to combined data without navigator 
%   @param  dst_navi         path to combined data with navigator correction 
    
%   Author: Martin Soellradl
%   Department of Neurology, Medical University of Graz, Graz, Austria
%   email:martin.soellradl@medunigraz.at
%   Website: http://www.neuroimaging.com
%   Januray 2020; Last revision: 02-Januray-2020

    if isunix; 
    	
    else
        warning('Unix with FSL Flirt required! Exit script'); 
    end

    if nargin < 5
       bloadExample = 0;  
    end

    %%  Load and correct data with the navigator echo 
    Navg = length(dat_path_c); 
    
    for j=1:Navg; 
        S_raw_path = [dst_raw, '/meas_', num2str(j), '/']; 
        S_corr_path =  [dst_navi, '/meas_', num2str(j), '/'];  
        [Scomb_corr(:,:,:,:,j), Scomb_raw(:,:,:,:,j)  ] = prepareAndCorrectRawData( dat_path_c{j}, S_raw_path, S_corr_path, bdetrend, nii_headers{j},bloadExample);
    end
    
    res = []; 
            
    %% unwrap phase data with prelude 
    
    dst_files{1} = dst_raw; 
    dst_files{2} = dst_navi;
    navi{1} = 'navi_off';  navi{2} = 'navi_on'; 
    
    for i=1:length(dst_files)
        for j=1:Navg; 
              unwrapped = [dst_files{i},'/meas_', num2str(j), '/phase_', navi{i}, '_unwrapped.nii.gz']; 
              if ~exist(unwrapped, 'file')
                  %path to mag and phase (4D)
             
                  phasevol = [dst_files{i},'/meas_', num2str(j), '/phase_', navi{i}, '.nii.gz'];  
                  absvol = [dst_files{i},'/meas_', num2str(j), '/mag_', navi{i}, '.nii.gz']; 

                  %destination unwrapped phase 
                  path = [dst_files{i},'/meas_', num2str(j), '/']; 


                  tmp = load_untouch_nii(phasevol); 
                  phase_tmp = double(tmp.img); 
                  %mask_path = [path_results, '/Mask/', 'mask_', file_id]; 

                  prelude_call = ['prelude ', '-p ', phasevol, ' -a ', absvol, ...
                      ' -o ',  unwrapped, ' -s']; 
                  system(prelude_call); 

                  %load unwrapped mask 
                  tmp = load_untouch_nii(unwrapped); 
                  phase_prelude{i,j} = double(tmp.img); 
              else 
                  tmp = load_untouch_nii(unwrapped); 
                  phase_prelude{i,j} = double(tmp.img); 
              end


        end
    end
    
    
    %% register images with FSL Flirt 

    for i=1:2
        
        %paths to raw data and navi corrected data
        if i==1
           path_data = [dst_raw]; 
           mag_name = 'mag_navi_off'; 
           phase_name = 'phase_navi_off_unwrapped'; 
        else
           path_data =  [dst_navi];
           mag_name = 'mag_navi_on'; 
           phase_name = 'phase_navi_on_unwrapped'; 
        end
        
        %register averages 
        ref_mag_path = [path_data, 'meas_1/', mag_name, '.nii.gz']; 
        ref_phase_path = [path_data, 'meas_1/', phase_name, '.nii.gz']; 
        for j=2:Navg; 
            in_path_mag =   [path_data, 'meas_', num2str(j), '/', mag_name, '.nii.gz'];  
            out_path_mag = [path_data, 'meas_', num2str(j), '/', mag_name, '_reg.nii.gz'];  

            in_path_phase = [path_data, 'meas_', num2str(j), '/', phase_name, '.nii.gz'];  
            out_path_phase = [path_data, 'meas_', num2str(j), '/', phase_name, '_reg.nii.gz'];  
            
            if ~exist(out_path_phase, 'file'); 
            
                out_Tmat= [path_data, 'meas_', num2str(j), '/' 'Tmat_', mag_name];

                options_flirt = ['-dof 6 -searchrx -15 15 -searchry -15 15 -searchrz -15 15  -coarsesearch 30 -finesearch 10 -cost mutualinfo ']; 
                flirt_call = ['flirt ', options_flirt,  '-in ', in_path_mag, ' -ref ', ref_mag_path,' -out ', out_path_mag, ' -omat ', out_Tmat]; %register mGRE images on T1w image 
                system(flirt_call);


                %transform phase of the second acquistion 
                flirt_call = ['flirt -in ', in_path_phase, ' -ref ', ref_mag_path, ...
                                    ' -applyxfm ', ' -init ', out_Tmat,' -out ', out_path_phase]; 
                system(flirt_call);            


                %transform mag again...apperently flirt just transforms first echo 
                flirt_call = ['flirt -in ', in_path_mag, ' -ref ', ref_mag_path, ...
                                    ' -applyxfm ', ' -init ', out_Tmat,' -out ', out_path_mag]; 
                system(flirt_call);  
               
            end

            %load data 
            tmp = load_untouch_nii(ref_mag_path); 
            Smag{i,1} = double(tmp.img); 
            tmp = load_untouch_nii(ref_phase_path); 
            phase{i,1} = double(tmp.img); 

            tmp = load_untouch_nii(out_path_mag); 
            Smag{i,j} = double(tmp.img); 
            tmp = load_untouch_nii(out_path_phase); 
            phase{i,j} = double(tmp.img); 
            
        end
        
        
    end
    
     phase1 = phase{2,1}; 
     phase2 = phase{2,2}; 
     
     phase_avg = (phase1 + phase2)/2;
     
    imagine(squeeze(phase1(:,:,:,2)), squeeze(phase2(:,:,:,2)), squeeze(phase_avg(:,:,:,2))); 
    
    
    res.Smag_sep = Smag; 
    res.phase_sep = phase;

    for i=1:2
        
       if i==1
           dst_avg = [dst_raw, '/meas_sep_avg/']; 
       else
           dst_avg = [dst_navi, '/meas_sep_avg/']; 
       end
        
       if ~exist(dst_avg, 'dir')
            mkdir(dst_avg); 
       end
     
       mag_avg_path = [dst_avg, '/mag_sep_avg.nii.gz']; 
       phase_avg_path = [dst_avg, '/phase_sep_avg.nii.gz']; 
        
       if ~exist(mag_avg_path, 'file'); 
           mag_avg_raw =  (Smag{i,1} + Smag{i,2})/2; 
           phase_avg_raw =  (phase{i,1} + phase{i,2})/2; 
           
           nii = nii_headers{i}; 
           nii.hdr.dime.datatype = 16; 
           nii.hdr.dime.bitpix = 32;  
           nii.img  = mag_avg_raw; 
           save_untouch_nii(nii, mag_avg_path); 
           
           
           nii.img  = phase_avg_raw; 
           save_untouch_nii(nii, phase_avg_path); 
        
       end
        
        if i==1
           tmp = load_untouch_nii(mag_avg_path); 
           res.mag_avg_raw = double(tmp.img); 
           tmp = load_untouch_nii(phase_avg_path); 
           res.phase_avg_raw = double(tmp.img); 
        else
           tmp = load_untouch_nii(mag_avg_path); 
           res.mag_avg_navi = double(tmp.img); 
           tmp = load_untouch_nii(phase_avg_path); 
           res.phase_avg_navi = double(tmp.img); 
        end


    
    end




end

