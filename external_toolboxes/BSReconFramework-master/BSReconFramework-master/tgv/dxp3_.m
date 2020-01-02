% x forward differentiation negative adjoint 
function  v= dxp3_(u)
    N = [size(u,1) 1 size(u,3)];
    v = cat(2, u(:,1:end-1,:), zeros(N)) ...
        - cat(2, zeros(N), u(:,1:end-1,:));
end