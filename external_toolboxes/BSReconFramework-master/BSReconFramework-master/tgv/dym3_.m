% y backwards differentiation negative adjoint 
function  v= dym3_(u)
    N = [1 size(u,2) size(u,3)];
    v = cat(1, zeros(N), -u(2:end,:,:)) ...
        + cat(1, u(2:end,:,:), zeros(N));
end