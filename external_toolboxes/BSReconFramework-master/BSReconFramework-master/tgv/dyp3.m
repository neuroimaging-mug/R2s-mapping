% y forward differentiation
function v = dyp3(u)
    N = [1 size(u,2) size(u,3)];
    v = cat(1, u(2:end,:,:) - u(1:end-1,:,:), zeros(N));
end

