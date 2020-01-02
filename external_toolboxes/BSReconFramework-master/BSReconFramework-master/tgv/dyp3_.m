% y forward differentiation negative adjoint 
function  v= dyp3_(u)
    N = [1 size(u,2) size(u,3)];
    v = cat(1, u(1:end-1,:,:), zeros(N)) ...
        - cat(1, zeros(N), u(1:end-1,:,:));
end
