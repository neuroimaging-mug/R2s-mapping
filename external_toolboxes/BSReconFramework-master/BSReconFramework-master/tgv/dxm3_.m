% x backwards differentiation negative adjoint 
function  v= dxm3_(u)
    N = [size(u,1) 1 size(u,3)];
    v = cat(2, zeros(N), -u(:,2:end,:)) ...
        + cat(2, u(:,2:end,:), zeros(N));
end
