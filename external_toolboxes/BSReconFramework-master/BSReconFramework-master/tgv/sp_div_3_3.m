function [y] = sp_div_3_3(x,Dx,Dy,Dz,dx,dy,dz)

[n,m,t,k] = size(x);
N = n*m*t;

y =  - reshape( ( Dx'*reshape(x(:,:,:,1),N,1) )./dx + ...
                ( Dy'*reshape(x(:,:,:,2),N,1) )./dy + ...
                ( Dz'*reshape(x(:,:,:,3),N,1) )./dz ,n,m,t );
		
