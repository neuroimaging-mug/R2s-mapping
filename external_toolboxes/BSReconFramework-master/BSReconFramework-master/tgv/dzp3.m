% z forward differentiation
function v = dzp3(u)
    N = [size(u,1) size(u,2) 1];
    v = cat(3, u(:,:,2:end) - u(:,:,1:end-1), zeros(N));
end
