function [ par_c ] = loadExampleSubjectInfoCell(  )
%LOADSUBJECTINFOCELL File containing file informations for each subject.

    %path to folder with all nii-files 
    par.path_pwd  = pwd; 

    %name of a specific acquistion. Here one with positive and negative
    %Gslice polarity. 
    par.acquisitions{1} = 'sinc_hanning_BWT_2d7_Tpulse_2ms_GsPos';
    par.acquisitions{2} = 'sinc_hanning_BWT_2d7_Tpulse_2ms_GsNeg';
    
    %Boolean if paramter a B1-map is available. In case not, calculations are
    %performed without B1. 
    par.bB1_map = 1;
    
    %if phase offset phi0z  is available the gradient Gphi0z can also be
    %considered (dw(t) = phi0z + gamma*Gz*z*t). 
    %With the here provided data it is not possible due to the
    %coil combination (phi0z cancels out)
    par.Gphi0z = 0; 
    
    %Here only caculations are performed for a single slice
    par.sgl_sag_slc = 0; 
    
    %Paramter for the brain extraction tool (BET). In case no BET/FSL is avaible
    %a simple mask by thresholding is generated.
    par.bbet = 1; 
    
    %Performs calculation only for the data with navigator if true
    par.navi_only = 1; 
    
    
    %detrends phase evolving of the navigator signal. 
    %This is due to the phase of the ADC. Phase increases linear with the
    %number of phase encoding steps and thus has to be removed to obtain
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

    B1_new_folder =[par.path_pwd, '/data_input/Example_2_R2s/', par.meas_id, '/B1_map/'];
    par.src_nii = [par.path_pwd, '/data_input/Example_2_R2s/', par.meas_id, '/nifti/'];
    par.src_dcm = [par.path_pwd, '/data_input/Example_2_R2s/', par.meas_id, '/dcm/'];
    src_raw_data = [par.path_pwd, '/data_input/Example_2_R2s/', par.meas_id, '/dat_file/'];
    
    
    %Boolean for the provided example data 
    par.bloadExample=1; 
    
    %Here useally the path to the SIEMENS data file is entered. Because of data 
    %protection instead of a dat file for the example a *.mat file is loaded 
    %with the coil data and a header with all necessary paramters. 
    par.dat_path{1} = [src_raw_data, 'GRE_pos_data_coil'];
    par.dat_path{2} = [src_raw_data, 'GRE_neg_data_coil'];
        
    par.nii_file{1} = [  par.src_nii , 'Gs_pos_myRF_normal_navi_18E_TR_2_5_s011.nii.gz'];
    par.nii_file{2} = [  par.src_nii , 'Gs_neg_myRF_normal_navi_18E_TR_2_5_s013.nii.gz'];
    
    par.bB1_map = 1;
    par.bB1_map_new = 1;
    par.B1_path = [B1_new_folder, 'B1Map_prospectivelySubsampled_subject_1_meas_MID00397_FID105235_BS_undersampled_R2s'];
    par_c{1} = par; 


end

