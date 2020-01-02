function BSrecon_all = h1_solve_cg(data,encop,mu,maxit,tol,imgdims, dz)

% TODO:  tune mu
%[ny,nx,ncoils] = size(data);
ny = imgdims(1);
nx = imgdims(2);
nz = imgdims(3);

A = laplace_matrix3DnewScale( ny,nx,nz, [],'sym', dz);

M = @(x) mu.*(encop'*(encop*x)) + A*x;
rhs = mu.*(encop'*data(:));

BSrecon = pcg(M,rhs,tol,maxit);
BSrecon_all = reshape(BSrecon,[ny,nx,nz]);

end

