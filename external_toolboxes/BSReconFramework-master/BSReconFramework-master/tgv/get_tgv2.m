

function [tgv] = get_tgv2(x,alph0,alph1,dx,dy,dz)
	
	%Get [ Du - v , E(v) ]
	x0 =  cat( 4 , fgrad_3( x(:,:,:,1), dx,dy,dz ) - x(:,:,:,2:4), sym_bgrad_3( x(:,:,:,2:4), dx,dy,dz ) );
	
    %Set norm
    tgv = alph1*(sum(sum(sum( norm_3( abs(x0(:,:,:,1:3)) ) ))) ) + alph0*(sum(sum(sum( norm_6( abs(x0(:,:,:,4:9)) ) ))) );

end
