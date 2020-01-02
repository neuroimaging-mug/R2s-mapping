function [ data ] = forward_mri3d(img, b1, mask)
% Forward operation: image to kspace 
[ny,nx,nz,ncoils] = size(b1);
data = zeros(ny,nx,nz,ncoils);

for coil=1:ncoils
    data(:,:,:,coil) = mask.*ifftn(img.*b1(:,:,:,coil)).*sqrt(ny*nx*nz);  

end

