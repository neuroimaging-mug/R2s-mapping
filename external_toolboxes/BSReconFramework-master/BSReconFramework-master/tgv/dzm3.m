% z backwards differentiation
function v = dzm3(u)
    N = [size(u,1) size(u,2) 1];
    v = cat(3, zeros(N), u(:,:,2:end) - u(:,:,1:end-1));
end