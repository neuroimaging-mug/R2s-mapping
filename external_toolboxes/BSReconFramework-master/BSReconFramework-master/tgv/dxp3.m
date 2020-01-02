function v = dxp3(u)
% x forward differentiation
    N = [size(u,1) 1 size(u,3)];
    v = cat(2, u(:,2:end,:) - u(:,1:end-1,:), zeros(N));
end

