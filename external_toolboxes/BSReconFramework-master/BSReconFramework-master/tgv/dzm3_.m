% z backwards differentiation negative adjoint 
function  v= dzm3_(u)
    N = [size(u,1) size(u,2) 1];
    v = cat(3, zeros(N), -u(:,:,2:end)) ...
        + cat(3, u(:,:,2:end), zeros(N));
end
