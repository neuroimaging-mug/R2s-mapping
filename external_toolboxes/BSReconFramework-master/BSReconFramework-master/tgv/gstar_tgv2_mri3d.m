function [g] = gstar_tgv2_mri3d(x, y, z, mri_obj, dx, dy, dz, ld)


	%F(Kx)(data, b1, mask
	g1 = abs( forward_mri3d ( x(:,:,:,1), mri_obj.b1, mri_obj.mask ) - mri_obj.data);	
    g1 = (ld/2)*sum(g1(:).^2);
	
	%F*(z)
	g2 = sum( mri_obj.data(:)'*z(:) ) + (1/(2*ld))*sum(abs(z(:)).^2);
	
	%G*(-Kx) (not for the real gap)
	g3 = - bdiv_3( y(:,:,:,1:3) , dx,dy,dz ) + reshape( backward_mri3d(z,mri_obj.b1,mri_obj.mask) , [ size(x,1), size(x,2), size(x,3)]);
	g3 = sum(abs(g3(:)));
	
	g4 = -y(:,:,:,1:3) - fdiv_3( y(:,:,:,4:9) , dx,dy,dz ) ;
	g4 = sum(abs(g4(:)));
	
	g =   g1 + g2 + g3 + g4;