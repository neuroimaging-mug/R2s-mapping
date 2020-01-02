
function [sig,tau] = steps_tgv2(x,sig,tau,K,dx,dy,dz)

	[n,m] = size(x);

	Kx = zeros(n,m,6);    
	
	%Get Kx
    	Kx = abs( cat( 4 ,  fgrad_3( x(:,:,:,1) , dx,dy,dz ) - x(:,:,:,2:4)	,...
    				sym_bgrad_3( x(:,:,:,2:4) ,  dx,dy,dz	) ) );
    				
 	Kx2 = abs( K( x(:,:,:,1)) );
 	
    	%Get |Kx|
    	nKx = sqrt(	sum(sum(sum( 	Kx(:,:,:,1).^2 + Kx(:,:,:,2).^2 + Kx(:,:,:,3).^2 + ...
    					Kx(:,:,:,4).^2 + Kx(:,:,:,5).^2 + Kx(:,:,:,6).^2 + ...
    					2*Kx(:,:,:,7).^2 + 2*Kx(:,:,:,8).^2 + 2*Kx(:,:,:,9).^2  ))) + ...
    			sum(		Kx2(:).^2 )   							);
    	
    	%Get |x|
    	nx = sqrt(sum(sum(sum(	abs(x(:,:,:,1)).^2 + abs(x(:,:,:,2)).^2 + abs(x(:,:,:,3)).^2 + abs(x(:,:,:,4)).^2  ))));
    	                    
    
    %Set |x| / |Kx|
    tmp = (nx/nKx);
    theta = 0.95;
    
    %Check convergence condition
    if sig*tau > tmp^2 %If stepsize is too large
        if theta^(2)*sig*tau < tmp^2 %Check if minimal decrease satisfies condition
            sig = theta*sig;
            tau = theta*tau;
        else                        %If not, decrease further
            sig = tmp;
            tau = tmp;
        end
    end
    
    %sig = tmp;
    %tau = tmp;
    
end
