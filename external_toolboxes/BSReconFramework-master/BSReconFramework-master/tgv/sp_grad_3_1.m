function [y] = sp_grad_3_1(x,Dx,Dy,Dz,dx,dy,dz)

[n,m,t] = size(x);

y = reshape( [ (Dx*x(:))/dx ; (Dy*x(:))/dy ; (Dz*x(:))/dz] , n,m,t,3);