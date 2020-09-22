
function [ res ] = performMonoexponentialFit( mag,mask,te,file_id, path_results,opts )
%PERFORMMONOEXPONENTIALFIT performs a monoexponential fit of standard
%gradient-echo data without z-shim
                 
                 
    % ----------------------------------------------------------------
    % Perform a monoexponential fit of the data 
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
    

end

