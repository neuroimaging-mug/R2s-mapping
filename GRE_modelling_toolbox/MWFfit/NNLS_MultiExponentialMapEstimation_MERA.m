function [MWF_map, M0_map,out1D] = NNLS_MultiExponentialMapEstimation_MERA(mag, te, mask, beta_reg)
%R2SMAPESTIMATION Estimates MWF estimation with 
%a non-negatvie least squares (NNLS) multi-compartment

%   Script that performs an R2s estimation from multi gradient echo data. 
%   @param      mag      Magnitude image [Nx, Ny, Nz, t]
%   @param      te       echo times in ms 
%   @param      mask    
%   @param      beta_reg regularization paramter for the NNLS fit     

% Author: Martin Soellradl
% Department of Neurology, Medical University of Graz, Graz, Austria
% email:martin.soellradl@medunigraz.at
% Website: http://www.neuroimaging.com
% Januray 2020; Last revision: 02-Januray-2020

    th_MWF = 0.025; %s


    if nargin < 4
        beta_reg = 0.001; 
    end
    
    
    fitting.regtyp = 'me';
    % % optional inputs for 'mg' fitting
    % fitting.widthgauss = 5;
    % fitting.numbergauss = 2;
    % fitting.T0gauss = [0.02 0.07]';
    % fitting.numberT0 = 5;

    % for conventional regularization, the regulization weighting must be
    % adjusted. This can be done manually (fitting.regadj = 'manual') or using
    % one of two automatic methods ('gcv' and 'erinc') -- see literature
    % references in code
    fitting.regadj = 'manual';
    fitting.regweight = beta_reg; % 0.000000001; % 0.001

    % graph the results or not. This input is irrelevant if
    % analysis.interactive = 'y'

    analysis.graph = 'n';

    % define the range of time constants to fit. Note that for 'mg' fitting,
    % this is the full range of the spectrum, but the lowest and highest mean
    % component time constants (echoed to display) cover a narrower domain
    fitting.rangeT = [2.5e-3 0.25];

    % define the number of time constants in the spectrum
    fitting.numberT = 200;

    % set the non-negative least square (NNLS) code. In most cases, the supplied
    % nnlsmex fuction is the fastest NNLS solver. You may need to compile this
    % for your system. Or you can specify the MATLAB fuction,
    % fitting.nnlscode = 'lsqnonneg'. 

    % You can automatically or manually extract a finite number of component
    % amplitudes and time constants. Test this out with the interactive GUI.
    analysis.extract = 'auto';
    analysis.numberextract = 2;


    te = te*1E-3; %must be in s
    [Nx, Ny, Nz, NE] = size(mag); 
    M0_map = zeros(Nx, Ny, Nz); 
    MWF_map = M0_map; 
    
    %prepare data for parfor
    data.t = te';
    for i=1:Nz
        IMGr = reshape(mag(:,:,i,:),Nx*Ny,NE);
        bwr = reshape(mask(:,:,i),Nx*Ny,1);
        X = IMGr(bwr > 0,:)'; % NOW X is an NE x Nvox matrix of decay data
        data.D = X;
        data_c{i} =  data;
    end

    parfor i= 1:Nz 
       disp(['slice i: ', num2str(i), ' out of', num2str(Nz)]); 
        
       [k] = find(mask(:,:,i) > 0); 

        data_tmp =  data_c{i}; 
        if ~isempty(data_tmp.D); 
            
            tic
                [out1D]=MERA(data_tmp,fitting,analysis);
            toc

            %The fitted signal amplitude at t = 0
            M0_vec = sum(out1D.S);

            %Estimate MWF 
            [~, idx_th] = min(abs( out1D.T - th_MWF)); 
            MWF_vec = sum(out1D.S(1:idx_th,:))./sum(out1D.S);


            MWF = zeros(Nx, Ny);   
            M0 = MWF; 
            for vx=1:length(k); 
                MWF(k(vx)) = MWF_vec(vx); 
                M0(k(vx)) =  M0_vec(vx); 
            end

            MWF_map(:,:,i) = MWF; 
            M0_map(:,:,i) = M0; 
        end
    end
end

