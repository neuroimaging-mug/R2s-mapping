% x backwards differentiation
function v = dxm3(u)
    N = [size(u,1) 1 size(u,3)];
    v = cat(2, zeros(N), u(:,2:end,:) - u(:,1:end-1,:));
end