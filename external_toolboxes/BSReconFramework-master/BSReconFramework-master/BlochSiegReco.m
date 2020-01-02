function [B1Map, flipAngleMap] = BlochSiegReco(par)


% Comments
% Reconstruction algorithm for B1+ maps from highly undersampled 
% Bloch-Siegert data based on the method described in ..............
% ............................................................ Zitat
% 
% inpur parameter: struct par ... containing input parameters
%           par.pulse      ... normalized shape of Bloch-Siegert pulse
%           par.Tpulse     ... duration of BS-pulse in ms
%           par.deltaf     ... ofresonance frequency of BS-pulse in Hz
%           par.deltaOmega ... 2*pi*par.deltaf
%           par.alphaBS    ... onresonant flip angle of BS-pulse in deg
%           par.gamma      ... gyromagnetic ratio 
%           par.deltat     ... dwell time of BS-pulse
%           par.Kbs        ... pulse constant in rad/G²
%           par.impl       ... for 'CPU': use CPU implementation
%                              for 'GPU': use GPU implementation
%           par.pattern    ... undersampling pattern
%           par.NC         ... number of channels
%           par.dimY       ... number of phaseEncoding Lines
%           par.dimX       ... number of sampling points in frequency
%                              encoding direction
%           par.dimSlice   ... number of encodings in 2nd phase encoding
%                              direction (partitions)
%           par.y          ... undersampled k-space data
%           par.PhaseEncDir .. phase encoding direction
%           par.dx         ... Resolution in x-direction
%           par.dy         ... Resolution in y-direction         
%           par.dz         ... Resolution in z-direction
%           par.traj       ... trajectory:  'cart': only Cartesian implemented
%           par.loadCoilSens . 
%           par.CoilSensPath .
%           par.CoilSensData .
%           par.lambda     ... Regularization parameter for TGV reconstruction
%           par.mu         ... Regularization parameter for H1 reconstruction
%
%
% output parameters:
%           B1Map          ... final 3D B1 map in µT
%           flipAngleMap   ... corresponding flip angle correction map for
%                              onresonant excitation -> 100% correspond to
%                              nominal flip anlge

%% normalize data and calculate/load coil sensitivity maps
[ny,nx,nz,~,ncoils] = size(par.y);

% used for data normalization
calibsizeXY = 16;
calibsizeZ  = 16;
filterKernel1DXY = hamming(calibsizeXY);
filterKernel1DZ = hamming(calibsizeZ);
filterKernel = repmat(filterKernel1DXY*filterKernel1DXY', [1,1,calibsizeZ]);
filterKernel = filterKernel .* repmat(reshape(filterKernel1DZ, [1,1,calibsizeZ]), [calibsizeXY, calibsizeXY, 1]);
fmask = zeros(ny,nx,nz);
fmask((ny-calibsizeXY)/2+1:(ny+calibsizeXY)/2, (nx-calibsizeXY)/2+1:(nx+calibsizeXY)/2, (nz-calibsizeZ)/2+1:(nz+calibsizeZ)/2) = filterKernel;
fmask = repmat(fmask,[1, 1, 1, ncoils]);
 
if strcmp(par.coilSens, 'calcFromFullData')
    disp ('Calculating CoilSensitivity maps')
    calibdata = ifft3c(par.CoilSensData);
    [~, smaps] = walsh_sens_3d_slice( calibdata );
elseif strcmp(par.coilSens, 'loadExternal')
    disp ('Loading CoilSensitivity maps')
    disp(par.CoilSensPath);
    load(par.CoilSensPath);
else
    error('No coil sensitivity maps defined!')
end


dataLowRes = squeeze(par.y(:,:,:,1,:));
lowres = sqrt(sum(abs(ifft3c(fmask.*dataLowRes)).^2,4));
dscale = 100/max(lowres(:));
par.y = par.y.*dscale;

disp('Initialization finished')


%% start reconstruction
if strcmp(par.impl, 'CPU')
    %% Perform calculation on CPU
    % 1) TGV recon data_1

    reduction = 1;
    alpha = 1;

    % define operator
    switch par.traj
        case 'cart'
            encop = ENC_OP(smaps,par.pattern(:,:,:,1),1,0);
        case 'radial'
            encop = ENC_OP(smaps,par.pattern,w,1);
    end
    
    K = @(x) encop*x;
    Kh = @(x) encop'*x;
    col = @(x) x(:);
    
    % test adjoint
    u = rand(ny,nx,nz);
    v = rand( sum(col(par.pattern(:,:,:,1)))*ncoils, 1 );
    disp(['Adjointness TGV Operator', num2str(dot(K(u),v(:)) - dot(u(:),Kh(v)))]);
    
    datain = (squeeze(par.y(:,:,:,1,:)));
    datain = datain(repmat(par.pattern(:,:,:,1),[1 1 1 ncoils]));
    
    lambda = 1/par.lambda*2.8824;
    
    TGV2_recon = tgv2solve_pd_3d_new(datain,encop, lowres(:),lambda, sqrt(3)*alpha, alpha, par.maxitTGV, reduction,[ny,nx,nz,ncoils], 1,1,par.dz);
    
    % 2) H1 recon data_2
    % define operator
    disp('TGV reconstruction finished! Starting pcg for H1 ...');
    
    switch par.traj
        case 'cart'
            encop2 = ENC_OP(repmat(TGV2_recon,[1 1 1 ncoils]).*smaps,par.pattern(:,:,:,2),1,0);
        case 'radial'
            encop2 = ENC_OP(repmat(TGV2_recon,[1 1 ncoils]).*smaps,par.pattern,w,1);
    end
    
    datain = (squeeze(par.y(:,:,:,2,:)));
    datain = datain(repmat(par.pattern(:,:,:,2),[1 1 1 ncoils]));
    
    phaseImage1 = h1_solve_cg(datain,encop2,par.mu,par.maxitH1,1e-6,[ny,nx,nz,ncoils], par.dz);
    disp('H1 finished');

elseif strcmp(par.impl, 'GPU')
    %% Perform calculation on GPU
    par_in.maxItTGV = par.maxitTGV;
    par_in.lambda   = par.lambda;
    par_in.alpha0   = sqrt(3);
    par_in.alpha1   = 1;
    par_in.dx       = par.dx;
    par_in.dy       = par.dy;
    par_in.dz       = par.dz;
    par_in.parfile  ='./gpu/default.cfg';
    par_in.maxItH1  = par.maxitH1;
    par_in.mu       = par.mu;
    par_in.relTol   = 1e-24;
    par_in.absTol   = 1e-6;
    
    [~, phaseImage1, ~ ] = recon_gpu(squeeze(par.y(:,:,:,1,:)), ...
        squeeze(par.y(:,:,:,2,:)), par_in, par.pattern(:,:,:,1), par.pattern(:,:,:,2), ...
        smaps, lowres);
    
else
    error('Implementation method not defined!')
end

%% Calculate B1+ map
phaseImage = angle(phaseImage1);
B1Map = sqrt(phaseImage / par.Kbs);
flipAngleMap = B1Map * par.gamma*sum(par.pulse)*par.deltat/max(par.pulse) *180/pi / par.alphaBS *100;

end
