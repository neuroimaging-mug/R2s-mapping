function [data,mask] = setup_data3d(data,mask)

[n,m,l,ncoils] = size(data);


%Get chop variable
[x,y,z] = meshgrid(1:m,1:n,1:l);
chop = (-1).^(x+y+z+1);
clear x y z

%Shift mask and data for matlab fft , recenter image by chop
mask = ifftshift( mask );


if ~isempty(data)
    for coil=1:ncoils
        data(:,:,:,coil) = ifftshift(chop.*data(:,:,:,coil));
    end
    
else
    data=0;
end
