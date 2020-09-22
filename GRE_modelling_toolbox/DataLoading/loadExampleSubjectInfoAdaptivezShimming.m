function [par_c ] = loadExampleSubjectInfoAdaptivezShimming(  )
%LOADPARAMETERCONFIG loads all relevant config paramters for the adapative
%z-shim measurements. 
%path to folder with all nii-files 
    
    %path to folder with all nii-files 
    par_zShim.path_pwd  = pwd; 

    %% define pulse

    %myRF normal 
    pulse.type = 'sinc-hanning'; 
    pulse.Tpulse = 2;           
    %ms 
    pulse.BWT = 2.7;                %-
    pulse.k_pulse = 33.1510;        %�T/m
    pulse.GsPolarity = 'positive';
    
    %% Common parameters 
   
    par_zShim.bB1_map = 1; %set if B1 is available 
    par_zShim.B1_path = [];
    par_zShim.bbet = 1; %performs FSL bet brain extraction
    par_zShim.sgl_sag_slc = 0; %flag for single slice (speed-up)
    par_zShim.navi_only = 0; %performs calculations only with the navigator corrected data
    par_zShim.bdetrend = 1; %removes a linear phase offset of the navigator phase due to the ADC
    

    % -------------------------------------------------------------------
    % Subject 1
    % -------------------------------------------------------------------
    
    
    par_zShim.meas_id = 'subject_1';
   
    par_zShim.src_nii = [par_zShim.path_pwd, '/data_input/Example_5_R2s_zshim/', par_zShim.meas_id, '/nifti/'];
    par_zShim.src_dcm = [par_zShim.path_pwd, '/data_input/Example_5_R2s_zshim/', par_zShim.meas_id, '/dcm/'];
    src_raw_data = [par_zShim.path_pwd, '/data_input/Example_5_R2s_zshim/', par_zShim.meas_id, '/dat_file/'];
    par_zShim.T1_file = '';

    
    % par_zShim pamters of standard acuistion without zShim 
    par_zShim.acquisitions{1} = 'zShim_off';
    par_zShim.pulse_c{1} = pulse;
    par_zShim.dat_path{1} = [src_raw_data, 'GRE_data_coil_zshim_off'];
    
    % Paramters of Nam pattern Gc = 220µT/m
    par_zShim.acquisitions{2} = 'zShim_Gc_220';
    par_zShim.pulse_c{2} = pulse;
    par_zShim.dat_path{2} = [src_raw_data, 'GRE_data_coil_Gc_constant'];

    
    %  Paramters of Gcsplit max 
    par_zShim.acquisitions{3} = 'zShim_Gc_split';
    par_zShim.pulse_c{3} = pulse;
    par_zShim.dat_path{3} = [src_raw_data, 'GRE_data_coil_Gc_split'];
    
    
    % index of the nifti file dcmHeader
    par_zShim.idx_nii = [7,8,10]; 

    % zShim properties
    
    par_zShim.zShim.pattern{1} = '+A_-A_0_+B_-B'; 
    par_zShim.zShim.Gc{1} = 0; %(zero because of standard mGRE data)
   
    par_zShim.zShim.pattern{2} = '+A_-A_0_+B_-B'; 
    par_zShim.zShim.Gc{2} = 220; %µT/m
    
    par_zShim.zShim.pattern{3} = 'Gc_split_max'; ; 
    par_zShim.zShim.Gc {3} = NaN; %µT/m

    %Echo number after which the z-shim moments are applied
    par_zShim.zShim.start_echo = 4; 
    par_zShim.zShim.Gz_comp_path = [par_zShim.path_pwd, '/data_input/Example_4_zshimTable_prescan/',  par_zShim.meas_id, '/zShimTable/zShimTable_subject_1.txt'];
    
    
    %B1-map
    par_zShim.B1_path = [par_zShim.path_pwd, '/data_input/Example_5_R2s_zshim/', par_zShim.meas_id, '/B1_map/subject_1_B1_map.mat'];
    par_c{1} = par_zShim;

end

