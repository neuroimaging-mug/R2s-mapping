function [ TGV_recon, phaseImage, pdgap ] = recon_gpu( dataTGV, dataH1, par_in, maskTGV, maskH1, b1_in, u0_in)

%function [ g2, pdgap ] = recon_gpu( data, par_in, mask, b1_in, u0_in)
% simple export and import to gpu from matlab

switch nargin
    case 0
        disp('Provide: TGVData, H1Data, input Parameters, TGVPattern, H1Pattern, Sensitivity Profiles')
        return
    case 1
        disp('Provide: H1Data, input Parameters, TGVPattern, H1Pattern, Sensitivity Profiles')
        return
    case 2
        disp('Provide: input Parameters, TGVPattern, H1Pattern, Sensitivity Profiles')
        return
    case 3
        disp('Provide: TGVPattern, H1Pattern, Sensitivity Profiles')
        return
    case 4
        disp('Provide: H1Pattern, Sensitivity Profiles')
        return
    case 5
        disp('Provide: Sensitivity Profiles')
        return
    case 6
        u0_in = zeros(size(dataTGV,1),size(dataTGV,2),size(dataTGV,3));
        disp('Start init')
    otherwise
        disp('Start init')
end
        



%Regularization parameters for coil construction
method = 'BS_RECON';
%Stepsize
sigma = 1/3;
tau = 1/3;

%---------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------

% set zero output
comp1 = 0; comp2 = 0; pdgap = 0; b1 = 0; u0 = 0;

[n,m,l,ncoils] = size(dataTGV);


% % rescale
[dataTGV,maskTGV] = setup_data3d(dataTGV,maskTGV);
[dataH1,maskH1] = setup_data3d(dataH1,maskH1);

% write data to binary file
writebin_vector(permute(dataTGV,[2 1 3 4]),...
    './data.bin');
writebin_vector(permute(maskTGV,[2 1 3]),...
    './mask.bin');
writebin_vector(permute(dataH1,[2 1 3 4]),...
    './dataH1.bin');
writebin_vector(permute(maskH1,[2 1 3]),...
    './maskH1.bin');


writebin_vector(permute(b1_in,[2 1 3 4]),...
    './b1.bin');
writebin_vector(permute(u0_in,[2 1 3]),...
    './u0.bin');


recon_cmd=['avionic -i ',num2str(par_in.maxItTGV),...
    ' --tgv2_3D.dx=',num2str(par_in.dx),...
    ' --tgv2_3D.dy=',num2str(par_in.dy),...
    ' --tgv2_3D.dz=',num2str(par_in.dz),...
    ' --tgv2_3D.lambda=',num2str(par_in.lambda),...
    ' --tgv2_3D.alpha0=',num2str(par_in.alpha0),...
    ' --tgv2_3D.alpha1=',num2str(par_in.alpha1),...
    ' --tgv2_3D.sigma=',num2str(sigma),...
    ' --tgv2_3D.tau=',num2str(tau),...
    ' --h1.relTol=', num2str(par_in.relTol),...
    ' --h1.absTol=', num2str(par_in.absTol),...
    ' --h1.mu=', num2str(par_in.mu),...
    ' --h1.maxIt=' num2str(par_in.maxItH1),...
    ' -m ',method,' -p ',par_in.parfile,...
    ' --dataBS ./dataH1.bin --maskBS ./maskH1.bin --finalOutBS ./resultH1.bin -d ', ...
    num2str(m),':',num2str(n),':',num2str(l),':',...
    num2str(m),':',num2str(n),':',num2str(l),':',num2str(ncoils), ...
    ':1  -u ./u0.bin -s ./b1.bin ./data.bin ./mask.bin ',...
    './resultTGV.bin'];


display(recon_cmd);

% run reconstruction
unix(recon_cmd);

% read results
TGV_recon = readbin_vector('./resultTGV.bin');
TGV_recon = permute(reshape(TGV_recon,[m,n,l]),[2 1 3]);

phaseImage = readbin_vector('./resultH1.bin');
phaseImage = permute(reshape(phaseImage,[m,n,l]),[2 1 3]);


if exist(['./PDGap'])==2
    pdgap = readbin_vector('./PDGap');
    pdgap = abs(pdgap);
    unix('rm ./PDGap');
end



% clean up
unix('rm ./u0.bin ./b1.bin ./resultTGV.bin ./resultH1.bin ./data.bin ./mask.bin ./dataH1.bin ./maskH1.bin');

end
