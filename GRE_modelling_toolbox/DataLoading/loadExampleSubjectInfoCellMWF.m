function [ par_c ] = loadExampleSubjectInfoCellMWF()
%LOADSUBJECTCONFIGMWF Loads paramters of all subjects for MWF estimation. 
    
%   Author: Martin Soellradl
%   Department of Neurology, Medical University of Graz, Graz, Austria
%   email:martin.soellradl@medunigraz.at
%   Website: http://www.neuroimaging.com
%   Januray 2020; Last revision: 02-Januray-2020


    
    par.path_pwd  = pwd; 
    
    %name of the acquistion
    par.acquisition = 'RF_fast_TR_2s_2avg';

    % unique id of the subject 
    par.meas_id =  'subject_1';
    
    % folder for B1-map 
    B1_folder = [par.path_pwd, '/data_input/Example_3_MWF/', par.meas_id, '/B1_map/'];

    par.src_nii = [par.path_pwd, '/data_input/Example_3_MWF/', par.meas_id, '/nifti/'];
    par.nii_file = [  par.src_nii, 'myRF_fast_navi_28E_Gs_pos_TR_2s_1__s003.nii.gz']; 
    par.src_dcm = [par.src_nii, 'dcm'];
    src_raw_data =  [par.path_pwd, '/data_input/Example_3_MWF/', par.meas_id, '/dat_file/']; 
   
    %Boolean for the provided example data 
    par.bloadExample=1; 
    
    %Here useally the path to the SIEMENS data file is entered. Instead of 
    % a dat file for the example a *.mat file is loaded with the coil data 
    %and a header with all necessary paramters because of data protection. 
    par.dat_avg_path{1} = [src_raw_data, 'GRE_MWF_data_coil'];
    par.dat_avg_path{2} = [src_raw_data, 'GRE_MWF_flipped_RO_data_coil'];
    
%     par.dat_avg_path{1} = [src_raw_data, 'meas_MID00391_FID105229_myRF_fast_navi_28E_Gs_pos_TR_2s_1_2x1_2x4mm_BW_500'];
%     par.dat_avg_path{2} = [src_raw_data, 'meas_MID00392_FID105230_myRF_fast_navi_flipped_RO_28E_Gs_pos_TR_2s_1_2x1_2x4mm_BW_500'];
    par.bB1_map = 1; 
    par.B1_new = 1; 
    par.B1_path = [B1_folder, '/B1Map_prospectivelySubsampled_subject_1_meas_MID00394_FID105232_BS_undersampled_MWF'];
    par.idx_nii = [2,5]; 
    
    
    %detrends phase evolving of the navigator signal. 
    %This is due to the phase of the ADC. Phase increases linear with the
    %number of phase encoding steps and thus has to be romved to get
    %phase flucuations for phase corretion. 
    par.bdetrend = 1; 
    

    par_c{1} = par; 


    


end

