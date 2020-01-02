function [ S_comb ] = coilCombinationMultiEchoLuo( raw_data )
%COILCOMBINATIONMULTIECHO_ Combines multi gradient echo data based on the
%the work fo Luo NeurImage 2012: 
% Gradient Ech%Plural Contrast ImagingSignal model and derived contrasts: 
% T2*, T1,Phase, SWI, T1f, FST2*and T2*-SWI
%
% Author: Martin Soellradl
% Department of Neurology, Medical University of Graz, Graz, Austria
% email:martin.soellradl@medunigraz.at
% Website: http://www.neuroimaging.com
% Januray 2020; Last revision: 02-Januray-2020

    Nx = size(raw_data,1); 
    Nphase = size(raw_data,2); 
    Nz = size(raw_data,3); 
    Nechoes = size(raw_data,4);
    Mcoils = size(raw_data,5); 

    for m=1:Mcoils
      disp(['Iff2c coil image m ', num2str(m), ' out of ', num2str(Mcoils)]); 
      for n=1:Nechoes; 
            coil_img(:,:,:,n,m) = ifft2c(squeeze(raw_data(:,:,:,n,m))); 
      end
    end


    noise_wdw = zeros(Nx, Nphase,Nz); 
    noise_wdw(1:10,1:10, :,:,:) = 1;
    noise_wdw(noise_wdw == 0) = NaN; 


    %Estiamte sigma for each coil
    for n=1:Nechoes 
        for m=1:Mcoils
            sigma(:,n,m) = squeeze(nanstd(nanstd(abs(coil_img(:,:,:,n, m)).*noise_wdw))); 
        end
    end

    %average sigmas of alle echoes
    sigma_avg = squeeze(mean(sigma,2));

    %sum up all chanels
    sigma_sum = sum(sigma_avg,2)./Mcoils;
    for m=1:Mcoils
       lambda(:,m) = sigma_sum./sigma(:,m); 
    end

    S_comb = zeros(Nx, Nphase, Nz, Nechoes); 
    
    %Combine images
    for slc=1:Nz; 
        disp(['Coil combination slice i ', num2str(slc), ' out of ', num2str(Nz)]); 
        for n=1:Nechoes
           sum_tmp = zeros(Nx, Nphase); 
           for m=1:Mcoils
             sum_tmp = sum_tmp + 1.*conj(coil_img(:,:,slc,1,m)).*coil_img(:,:,slc,n,m); %lambda(slc,m)
           end
           S_comb(:,:,slc,n) = sum_tmp./Mcoils; 
        end
    end



end

