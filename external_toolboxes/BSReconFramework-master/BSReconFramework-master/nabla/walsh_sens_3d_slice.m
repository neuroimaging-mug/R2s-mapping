function [recon,cmap, maxcoil1] = walsh_sens_3d_slice(yn,rn,norm, maxcoil)

% Reconstruction of array data and computation of coil sensitivities based
% on: a) Adaptive Reconstruction of MRI array data, Walsh et al. Magn Reson
% Med. 2000; 43(5):682-90 and b) Griswold et al. ISMRM 2002: 2410
%
%	Input:
%	yn: array data to be combined [ny, nx, nz, nc].
%	rn: data covariance matrix [nc, nc].
%	norm: =1, normalize image intensity
%
%	Output:
%	cmap: estimated coil sensitivity maps [ny, nx, nz, nc].

yn=permute(yn,[4,1,2,3]);
[nc,ny,nx,nz]=size(yn);
if nargin<3, norm=0; end
if nargin<2, rn=eye(nc);end

% find coil with maximum intensity for phase correction
if nargin<4
    [~,maxcoil]=max(sum(sum(sum(permute(abs(yn),[4 3 2 1])))));
end


bs1=8;  %x-block size
bs2=8;  %y-block size
st=1;   %increase to set interpolation step size

recon=zeros(ny,nx);

parfor z=1:nz
    yn_temp = yn(:,:,:,z);
    wsmall=zeros(nc,floor(ny./st),floor(nx./st));
    cmapsmall=zeros(nc,floor(ny./st),floor(nx./st));
    
    for x=st:st:nx
        for y=st:st:ny
            %Collect block for calculation of blockwise values
            ymin1=max([y-bs1./2 1]);
            xmin1=max([x-bs2./2 1]);
            % Cropping edges
            ymax1=min([y+bs1./2 ny]);
            xmax1=min([x+bs2./2 nx]);
            
            ly1=length(ymin1:ymax1);
            lx1=length(xmin1:xmax1);
            m1=reshape(yn_temp(:,ymin1:ymax1,xmin1:xmax1),nc,lx1*ly1);
            
            m=m1*m1'; %signal covariance
            
            % eignevector with max eigenvalue for optimal combination
            [e,v]=eig(inv(rn)*m);
            
            v=diag(v);
            [~,ind]=max(v);
            
            mf=e(:,ind);
            mf=mf/(mf'*inv(rn)*mf);
            normmf=e(:,ind);
            
            % Phase correction based on coil with max intensity
            %mf=mf.*exp(-j*angle(mf(maxcoil)));
            %normmf=normmf.*exp(-j*angle(normmf(maxcoil)));
            
            wsmall(:,y./st,x./st)=mf;
            cmapsmall(:,y./st,x./st)=normmf;
        end
    end
    
    % Interpolation of weights upto the full resolution
    % Done separately for magnitude and phase in order to avoid 0 magnitude
    % pixels between +1 and -1 pixels.
    for i=1:nc
        wfull(i,:,:,z)=conj(imresize(squeeze(abs(wsmall(i,:,:))),[ny nx],'bilinear').*exp(j.*imresize(angle(squeeze(wsmall(i,:,:))),[ny nx],'nearest')));
        cmap(i,:,:,z)=imresize(squeeze(abs(cmapsmall(i,:,:))),[ny nx],'bilinear').*exp(j.*imresize(squeeze(angle(cmapsmall(i,:,:))),[ny nx],'nearest'));
        %wfull(i,:,:)=conj(imageresize(squeeze(abs(wsmall(i,:,:))),[ny nx],1).*exp(j.*imageresize(angle(squeeze(wsmall(i,:,:))),[ny nx],0)));
        %cmap(i,:,:)=imageresize(squeeze(abs(cmapsmall(i,:,:))),[ny nx],1).*exp(j.*imageresize(squeeze(angle(cmapsmall(i,:,:))),[ny nx],0));
        
    end
end
% Phase correction based on coil with max intensity
cmap=cmap.*repmat(exp(-1i*angle(cmap(maxcoil,:,:,:))),[nc 1 1 1]);

recon=squeeze(sum(wfull.*yn));   %Combine coil signals.
% normalization proposed in the abstract by Griswold et al.
if norm, recon=recon.*squeeze(sum(abs(cmap))).^2; end

cmap=permute(cmap,[2,3,4,1]);
maxcoil1 = maxcoil;
