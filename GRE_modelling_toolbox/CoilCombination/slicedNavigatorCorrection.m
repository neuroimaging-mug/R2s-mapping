function [ S_corr, S_raw,  diff_angle, fit_res] = slicedNavigatorCorrection( raw_data, te, bdetrend)
%SLICEDNAVIGATORCORRECTION Pipeline for correction of fo fluctuations with
%a navigator echo. 
%Function allows a slice-wise correciton to save memory..
%   @param         raw_data        k- space raw data [Nread, Nphase, Nslices, Ncoils]. 
%   @param         te              Echo times [abritrary].
%   @param         bdetrend        Optional detrend paramater. In the old
%   navigator sequence version the ADC phase is linear increased. This issue
%   was resolved in the new one by setting the ADC phase to 0. 
%
% Author: Martin Soellradl
% Department of Neurology, Medical University of Graz, Graz, Austria
% email:martin.soellradl@medunigraz.at
% Website: http://www.neuroimaging.com
% Januray 2020; Last revision: 02-Januray-2020

    if nargin < 3
       bdetrend = 1;  
    end
    
    %Extract navigator data (last echo)
    raw_navi = raw_data(:,:,:,end,:);


    Nread = size(raw_data,1); 
    Nphase = size(raw_data,2); 
    Nslice = size(raw_data,3); 
    Nechoes = size(raw_data,4);
    Ncoils = size(raw_data,5); 

    %Estimate projection of each k-line 
    navi = fftshift(ifft(ifftshift(raw_navi,1), [], 1),1);
    
    clear raw_navi; 
    
    %create a mask 
    mask = nan(size(navi)); 
    mask(abs(navi) >= max(abs(navi(:))).*0) = 1; %mask does not further improve results->skipp it. 


    %subtract first navigator echo from the others by complex division 
    % navi_diff = navi_phase./repmat(navi_mask(:,1,:,:), [1, size(navi_mask,2),1,1]); 


    %try to unrwap navigator phase and then do the subtraction
    navi_angle = angle(navi); 
    
    navi_diff = navi.*repmat(conj(navi(:,1,:,:,:)), [1, Nphase, 1,1,1]); 
    navi_diff = navi_diff.*mask; 
    
    clear navi_angle 


    mean_navi = angle(nanmean(navi_diff.*mask,1)); 
    mean_navi_uwrp = unwrap(mean_navi,[], 2); 
    

    mean_navi(isnan(mean_navi)) = 0; 
    
    
    %do some stat
    sq_navi_unwrp = double(squeeze(mean_navi_uwrp)); 
    fit_res = zeros(2, Ncoils); 
    nx = 1:Nphase; 
    for coil=1:Ncoils
           ft = fittype( 'poly1' );
           phi_vec = sq_navi_unwrp(:,coil); 
           phi_fit = phi_vec(~isnan(phi_vec)); 
           nx_fit = nx(~isnan(phi_vec)); 

           if length(nx_fit) > 2
               % Fit model to data.
               [fitresult, gof] = fit(nx_fit(:), phi_fit(:), ft );
                fit_res(:, coil) = [fitresult.p1, fitresult.p2]'; 
           else
           end
    end
    
    
    
    mean_navi_uwrp(isnan(mean_navi_uwrp)) = 0; 
    
    if bdetrend == 1
        navi_det = zeros(size(mean_navi)); 
        for i=1:size(mean_navi,5)
            navi_det(:,:,1,1,i) = permute(detrend(permute(mean_navi_uwrp(:,:,1,1,i),[2,1,3,4,5])),[2,1,3,4,5]); 
        end
        diff_angle = repmat(navi_det,[Nread, 1, 1, 1, 1]); 
    else
        diff_angle = repmat(mean_navi_uwrp, [Nread, 1,1,1]);
    end
 
    %Correct k-space data with navi diff
    raw_corr = zeros(Nread, Nphase, Nslice, Nechoes-1, Ncoils);
    for i = 1:length(te)-1
        raw_corr(:,:,:,i,:) = raw_data(:,:,:,i,:).*exp(-1j*diff_angle.*te(i)/te(end));
    end

    [ S_corr ] = coilCombinationMultiEchoLuo( raw_corr );    
    [ S_raw ] = coilCombinationMultiEchoLuo( raw_data(:,:,:,1:end-1, :));
    %reorder data to match the nfiti files
    S_corr = imrotate(flip(S_corr,2),-90); 
    S_raw = imrotate(flip(S_raw,2),-90); 
    
end

