% z forward differentiation negative adjoint 
function  v= dzp3_(u)
    N = [size(u,1) size(u,2) 1];
    v = cat(3, u(:,:,1:end-1), zeros(N)) ...
        - cat(3, zeros(N), u(:,:,1:end-1));
end