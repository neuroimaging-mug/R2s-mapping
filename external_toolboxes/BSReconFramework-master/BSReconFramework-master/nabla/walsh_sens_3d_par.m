function [recon,cmap]=walsh_sens_3d_par(yn,rn,norm)

% Reconstruction of array data and computation of coil sensitivities based 
% on: a) Adaptive Reconstruction of MRI array data, Walsh et al. Magn Reson
% Med. 2000; 43(5):682-90 and b) Griswold et al. ISMRM 2002: 2410
%-------------------------------------------------------------------------
%	Input:
%	yn: array data to be combined [nz,ny, nx, nc]. 
%	rn: data covariance matrix [nc, nc].
%	norm: =1, normalize image intensity
%
%	Output:
%	recon: reconstructed image [nz,ny, nx].
%	cmap: estimated coil sensitivity maps [nz,ny, nx, nc].
%--------------------------------------------------------------------------
% Ricardo Otazo
% CBI, New York University
%--------------------------------------------------------------------------
%

yn=permute(yn,[4,1,2,3]);
[nc,nz,ny,nx]=size(yn);
if nargin<3, norm=0; end
if nargin<2, rn=eye(nc);end

% sum(sum(sum(permute(abs(yn),[4 3 2 1]))))
% find coil with maximum intensity for phase correction
[~,maxcoil]=max(sum(sum(sum(abs(yn),4),3),2));   

bs1=8;  %x-block size
bs2=8;  %y-block size
bs3=8;  %z-block size
st=1;   %increase to set interpolation step size

wsmall=zeros(nc,round(nz),round(ny),round(nx));
cmapsmall=zeros(nc,round(nz),round(ny),round(nx));


parfor x=1:nx
    if mod(x,10)==1
    disp(['x: ', num2str(x)])
    end
for y=1:ny
for z=1:nz    
    %Collect block for calculation of blockwise values
    ymin1=max([y-bs1./2 1]);                   
    xmin1=max([x-bs2./2 1]);                  
    zmin1=max([z-bs3./2 1]);                  
    % Cropping edges
    ymax1=min([y+bs1./2 ny]);                 
    xmax1=min([x+bs2./2 nx]);                  
    zmax1=min([z+bs3./2 nz]);
    
    ly1=length(ymin1:ymax1);
    lx1=length(xmin1:xmax1);
    lz1=length(zmin1:zmax1);
    m1=reshape(yn(:,zmin1:zmax1,ymin1:ymax1,xmin1:xmax1),nc,lz1*lx1*ly1);
      
    m=m1*m1'; %signal covariance
      
    % eignevector with max eigenvalue for optimal combination
    [e,v]=eig(inv(rn)*m);                    
                                               
    v=diag(v);
    [mv,ind]=max(v);
      
    mf=e(:,ind);                      
    mf=mf/(mf'*inv(rn)*mf);               
    normmf=e(:,ind);
    
    % Phase correction based on coil with max intensity
    mf=mf.*exp(-j*angle(mf(maxcoil)));        
    normmf=normmf.*exp(-j*angle(normmf(maxcoil)));

    wsmall(:,z,y,x)=mf;
    cmapsmall(:,z,y,x)=normmf;
end
end
end
wsmall = single(wsmall);
cmapsmall = single(cmapsmall);
recon=zeros(nz,ny,nx);

% Interpolation of weights upto the full resolution
% Done separately for magnitude and phase in order to avoid 0 magnitude 
% pixels between +1 and -1 pixels.
% if matlabpool('size') ~=0
%     parfor i=1:nc
%             wfull(i,:,:,:)=conj(resize_M(squeeze(abs(wsmall(i,:,:,:))),[nz ny nx],'bilinear').*exp(j.*resize_M(angle(squeeze(wsmall(i,:,:,:))),[nz ny nx],'nearest')));
%             cmap(i,:,:,:)=resize_M(squeeze(abs(cmapsmall(i,:,:,:))),[nz ny nx],'bilinear').*exp(j.*resize_M(squeeze(angle(cmapsmall(i,:,:,:))),[nz ny nx],'nearest'));
%     end
% else
    for i=1:nc
            wfull(i,:,:,:)=conj(resize_M(squeeze(abs(wsmall(i,:,:,:))),[nz ny nx],'bilinear').*exp(j.*resize_M(angle(squeeze(wsmall(i,:,:,:))),[nz ny nx],'nearest')));
            cmap(i,:,:,:)=resize_M(squeeze(abs(cmapsmall(i,:,:,:))),[nz ny nx],'bilinear').*exp(j.*resize_M(squeeze(angle(cmapsmall(i,:,:,:))),[nz ny nx],'nearest'));
    end
% % end
cmap = single(cmap);
wfull=single(wfull);
recon=squeeze(sum(wfull.*yn));   %Combine coil signals. 
% normalization proposed in the abstract by Griswold et al.
if norm, recon=recon.*squeeze(sum(abs(cmap))).^2; end

cmap=permute(cmap,[2,3,4,1]);
