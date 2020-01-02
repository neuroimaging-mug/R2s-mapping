function [ g2, pdgap ] = ...
    tgv3d_matlab_gpu( data, par_in, mask, b1_in, u0_in)
% simple export and import to gpu from matlab

if nargin < 4
    b1_in = [];
end

if nargin < 3
    mask = squeeze(data(:,:,1,:))~=0;
end

%Regularization parameters for coil construction
method = 'TGV2_3D';
stop_par = 500;
parfile = '/home/dieheart/workspace/AVIONIC/CUDA/config/default.cfg';

lambda = 1;
alpha0 = sqrt(3);
alpha1 = 1;
dx = 1;
dy = 1;
dz = 1;


%Stepsize
sigma = 1/3;
tau = 1/3;

%Read parameter-------------------------------------------------------------------------
%Input: par_in--------------------------------------------------------------------------
%Generate list of parameters
vars = whos;
for l=1:size(vars,1)
    par_list{l,1} = vars(l).name;
end
%Set parameters according to list
for l=1:size(par_in,1);
    valid = false;
    for j=1:size(par_list,1); if strcmp(par_in{l,1},par_list{j,1})
            valid = true;
            eval([par_in{l,1},'=','par_in{l,2}',';']);
        end; end
    if valid == false; warning(['Unexpected parameter at ',num2str(l)]); end
end
%---------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------

% set zero output
comp1 = 0; comp2 = 0; pdgap = 0; b1 = 0; u0 = 0;

[n,m,l,ncoils] = size(data);


% % rescale
[data,mask] = setup_data3d(data,mask);
 
b1_in = image_shift3d(b1_in);
u0_in = image_shift3d(u0_in);

% write data to binary file
writebin_vector(permute(data,[2 1 3 4]),...
    ['./data.bin']);
writebin_vector(permute(mask,[2 1 3]),...
    ['./mask.bin']);


writebin_vector(permute(b1_in,[2 1 3 4]),...
    ['./b1.bin']);
writebin_vector(permute(u0_in,[2 1 3]),...
    ['./u0.bin']);


recon_cmd=['/home/dieheart/workspace/AVIONIC/CUDA/bin/avionic -i ',num2str(stop_par),...
    ' --tgv2_3D.dx=',num2str(dx),...
    ' --tgv2_3D.dy=',num2str(dy),...
    ' --tgv2_3D.dz=',num2str(dz),...
    ' --tgv2_3D.lambda=',num2str(lambda),...
    ' --tgv2_3D.alpha0=',num2str(alpha0),...
    ' --tgv2_3D.alpha1=',num2str(alpha1),...
    ' --tgv2_3D.sigma=',num2str(sigma),...
    ' --tgv2_3D.tau=',num2str(tau),...
    ' -m ',method,' -p ',parfile,...
    ' -d ', ...
    num2str(m),':',num2str(n),':',num2str(l),':',...
    num2str(m),':',num2str(n),':',num2str(l),':',num2str(ncoils), ...
    ':1  -u ./u0.bin -s ./b1.bin ./data.bin ./mask.bin ./result.bin'];


display(recon_cmd);

% run reconstruction
unix(recon_cmd);

% read results
g2 = readbin_vector('./result.bin');
g2 = permute(reshape(g2,[m,n,l]),[2 1 3]);

%g2 = g2./dscale;

g2 = image_shift3d(g2);

if exist(['./PDGap'])==2
    pdgap = readbin_vector('./PDGap');
    pdgap = abs(pdgap);
    unix('rm ./PDGap');
end



% clean up
unix('rm ./u0.bin ./b1.bin ./result.bin ./data.bin ./mask.bin');

end
