function [ img ] = backward_mri3d(data, b1, mask)
% Backward operation: kspace to image 
[ny,nx,nz,ncoils] = size(b1);
img = zeros(ny,nx,nz);

for coil=1:ncoils
    img = img + conj(b1(:,:,:,coil)).*fftn(data(:,:,:,coil).*mask)./sqrt(ny*nx*nz);    
end

