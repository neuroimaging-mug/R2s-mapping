function [ par_c ] = loadExampleSubjectInfoCell(  )
%LOADSUBJECTINFOCELL File containing file informations for each subject.

    %path to folder with all nii-files 
    par.path_pwd  = pwd; 

    %name of a specific acquistion. Here one with positive and negative
    %Gslice polarity. 
    par.acquisitions{1} = 'sinc_hanning_BWT_2d7_Tpulse_2ms_GsPos';
    par.acquisitions{2} = 'sinc_hanning_BWT_2d7_Tpulse_2ms_GsNeg';
    
    %Paramter if B1-map is available. In case not, calculations are
    %performed without B1. 
    par.bB1_map = 1;
    
    %if phase offset phi0z is available the gradient Gphi0z can also be
    %considered. With the here provided data it is not possible due to the
    %coild combination (phi0z cancels out)
    par.Gphi0z = 0; 
    
    %Here only caculations are performed for a single slice
    par.sgl_sag_slc = 0; 
    
    %Paramter for the brain extraction tool (BET)
    par.bbet = 1; 
    
    %Performs calculation only for the data with navigator if true
    par.navi_only = 1; 
    
    %Parmater for siennax semgenation
    par.bSienna  = 0; 
    
    %detrends phase evolving of the navigator signal. 
    %This is due to the phase of the ADC. Phase increases linear with the
    %number of phase encoding steps and thus has to be romved to get
    %phase flucuations for phase corretion. 
    par.bdetrend = 1; 
    
    %corresponding T1-file (nifti) if available
    par.T1_file = '';
    
    %pulse properties of the RF pulse
    pulse.type = 'sinc-hanning'; 
    pulse.Tpulse = 2;           
    pulse.BWT = 2.7;                %-
    pulse.k_pulse = 33.1510;        %mm*µT/m
    pulse.GsPolarity = 'positive'; 
    pulse_c{1} = pulse; 

    %Polarity is chnaged for the second pulse
    pulse.GsPolarity = 'negative'; 
    pulse_c{2} = pulse;

    %store pulse cell array 
    par.pulse_c = pulse_c; 

  
    %% Subject specific input

    par.meas_id = 'subject_1';

    B1_new_folder =[par.path_pwd, '/data_input/Example_1_R2s/', par.meas_id, '/B1_map/'];
    par.src_nii = [par.path_pwd, '/data_input/Example_1_R2s/', par.meas_id, '/nifti/'];
    par.src_dcm = [par.path_pwd, '/data_input/Example_1_R2s/', par.meas_id, '/dcm/'];
    src_raw_data = [par.path_pwd, '/data_input/Example_1_R2s/', par.meas_id, '/dat_file/'];
    par.dat_path{1} = [src_raw_data, 'meas_MID00395_FID105233_Gs_pos_myRF_normal_navi_18E_TR_2_5s_1x1x3mm_BW_500'];
    par.dat_path{2} = [src_raw_data, 'meas_MID00396_FID105234_Gs_neg_myRF_normal_navi_18E_TR_2_5s_1x1x3mm_BW_500'];
        
    par.nii_file{1} = [  par.src_nii , 'Gs_pos_myRF_normal_navi_18E_TR_2_5_s011.nii.gz'];
    par.nii_file{2} = [  par.src_nii , 'Gs_neg_myRF_normal_navi_18E_TR_2_5_s013.nii.gz'];
    
    
    par.bB1_map = 1;
    par.bB1_map_new = 1;
    par.B1_path = [B1_new_folder, 'B1Map_prospectivelySubsampled_subject_14_meas_MID00397_FID105235_BS_undersampled_R2s'];
    par_c{1} = par; 


end

