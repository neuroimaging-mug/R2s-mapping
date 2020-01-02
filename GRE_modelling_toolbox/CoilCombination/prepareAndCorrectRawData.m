function [Scomb_corr, Scomb_raw  ] = prepareAndCorrectRawData( dat_path, uncorr_folder_path, corr_folder_path, bdetrend, nii_headers)
%PREPAREANDCORRECTRAWDATA Loads data-file and performs f0 correction of hte
%data. 
%   @param dat_path           path to the Siemens data file 
%   @param S_uncorr_path      path to the uncorrected combined data 
%   @param S_corr_path        path to the f0 corrected combined data 
%   @param bdetrend           boolian, detrends the phase of the navigator 
%                             signal if true 
%   @param nii_headers        nfti header where the combined data is store.
%                             If not provided a nifti is created with make_nii
%
% Author: Martin Soellradl
% Department of Neurology, Medical University of Graz, Graz, Austria
% email:martin.soellradl@medunigraz.at
% Website: http://www.neuroimaging.com
% Januray 2020; Last revision: 02-Januray-2020

    if ~exist(uncorr_folder_path, 'dir')
        mkdir(uncorr_folder_path); 
    end
     
    if ~exist(corr_folder_path, 'dir')
        mkdir(corr_folder_path); 
    end
    
    if nargin == 5
       bnii_hdr = 1;  
    else
       bnii_hdr = 0;   
    end

   
    if  exist( [uncorr_folder_path, 'mag_navi_off.nii.gz'], 'file') == 2
        tmp = load_untouch_nii([corr_folder_path, 'mag_navi_on.nii.gz']); 
        mag_corr = tmp.img; 
        tmp =  load_untouch_nii([corr_folder_path, 'phase_navi_on.nii.gz']); 
        phase_corr = tmp.img; 
        
        Scomb_corr = mag_corr.*exp(1i*phase_corr); 
        
        tmp = load_untouch_nii([uncorr_folder_path, 'mag_navi_off.nii.gz']); 
        mag_uncorr = tmp.img; 
        tmp =  load_untouch_nii([uncorr_folder_path, 'phase_navi_off.nii.gz']); 
        phase_uncorr = tmp.img; 
        Scomb_raw = mag_uncorr.*exp(1i*phase_uncorr); 
        
    else

        [data_coil, hdr, PhaseEncDir] = readRAW_BS_3D(dat_path);
        data_coil = single(data_coil); 
        Nslices = size(data_coil,3);
        for i=1:Nslices; 
            disp(['Slice is ', num2str(i), ' out of ', num2str(size(data_coil,3))]); 

            raw_data = (data_coil(:,:,i,:,:));
            
            nii = nii_headers; 
            Ny_nii = size(nii.img,1);
            Nphase = size(raw_data,2); 
            if Nphase ~= Ny_nii
                raw_data(:,Nphase+1:Ny_nii,:,:,:) = 0; 
            end


            %get echo times
            Nte = hdr.Phoenix.lContrasts; 
            te = cell2mat(hdr.Phoenix.alTE(1:Nte))/1000; 

            %get echo times
            Nte = hdr.Phoenix.lContrasts; 
            te = cell2mat(hdr.Phoenix.alTE(1:Nte))/1000; 

            [ S_corr(:,:,i,:), S_raw(:,:,i,:),  diff_angle_fit{i}] = slicedNavigatorCorrection( raw_data, te, bdetrend);

        end

%         if strcmp(PhaseEncDir, 'COL')
%             S_corr = flip(permute(S_corr, [2,1,3,4]),2);
%             S_raw = flip(permute(S_raw, [2,1,3,4]),2);
%             
%         end
        %% Put data in correct slice order 

        Scomb_corr = zeros(size(S_corr));
        Scomb_raw = zeros(size(S_raw));

        slc_order =  hdr.Config.chronSliceIndices; 

        for i=1:Nslices
            Scomb_corr(:,:,slc_order(i)+1,:) = S_corr(:,:,i, :); 
            Scomb_raw(:,:,slc_order(i)+1,:) = S_raw(:,:,i, :); 
        end


        %% Store data 
        
        Scomb_corr  =Scomb_corr*1E9; 
        Scomb_raw  = Scomb_raw*1E9; 
        
        %if header is avaiable us eheader
        if bnii_hdr == 1
            
            nii = nii_headers; 
            nii.img  = abs(Scomb_raw); 
            save_untouch_nii(nii, [uncorr_folder_path, 'mag_navi_off.nii.gz']); 
            nii.img  = angle(Scomb_raw); 
            save_untouch_nii(nii, [uncorr_folder_path, 'phase_navi_off.nii.gz']); 
                      
            nii.img  = abs(Scomb_corr);
            save_untouch_nii(nii, [corr_folder_path, 'mag_navi_on.nii.gz']); 
            nii.img  = angle(Scomb_corr); 
            save_untouch_nii(nii, [corr_folder_path, 'phase_navi_on.nii.gz']); 
            
        else


            tmp = make_nii(abs(Scomb_raw)); 

            %store header
            tmp.hdr.dime.datatype = 16; 
            tmp.hdr.dime.bitpix = 32;   

            save_nii(tmp, [uncorr_folder_path, 'mag_navi_off.nii.gz']); 

            tmp.img = angle(Scomb_raw); 
            save_nii(tmp, [uncorr_folder_path, 'phase_navi_off.nii.gz']); 


            tmp.img =  abs(Scomb_corr); 
            save_nii(tmp, [corr_folder_path, 'mag_navi_on.nii.gz']); 
            tmp.img = angle(Scomb_corr); 
            save_nii(tmp, [corr_folder_path, 'phase_navi_on.nii.gz']); 
        %     
        end
   end
end

