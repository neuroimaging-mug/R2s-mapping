% Computes the matrix A for the surface fitting problem (3)
% [nabla, valx, valy] = make_nabla(M,N,nz,mask)
% 
% Computes A for input images of size MxN
% nz is the z compenent of the normal map
% mask is a binary mask of foreground pixels
%
% Output:
% Matrix A for the surface fitting problem 
% valx ... vector of size M*N, 1 if unmasked depth measurement, 0 else
% valy ... vector of size 2*M*N, 1 if unmasked linear equation, 0 else
%
% i.e. A*z(valx) = b(valy)
%

function [nabla_, valx, valy] = make_nabla(M,N,mask)

row  = zeros(1,M*N*2);
col = zeros(1,M*N*2);
val  = zeros(1,M*N*2);

cnt = 1;

for y=1:M
  for x=1:N
      if mask(y,x) == 0
          continue;
      end
      if x < N && mask(x+1) == 1
          row(cnt) = y+(x-1)*M;
          col(cnt) = y+(x-1)*M;
          val(cnt) = -1;
          cnt = cnt+1;
    
          row(cnt) = y+(x-1)*M;
          col(cnt) = y+(x)*M;
          val(cnt) = 1;
          cnt = cnt+1;
      else
          row(cnt) = y+(x-1)*M;
          col(cnt) = y+(x-1)*M;
          val(cnt) = 1;
          cnt = cnt+1;
    
          row(cnt) = y+(x-1)*M;
          col(cnt) = y+(x-2)*M;
          val(cnt) = -1;
          cnt = cnt+1;
      end
  end
end
row = row(1:cnt-1);
col = col(1:cnt-1);
val = val(1:cnt-1);

Kxf = sparse(row,col,val,M*N,M*N);

row  = zeros(1,M*N*2);
col = zeros(1,M*N*2);
val  = zeros(1,M*N*2);

cnt = 1;
for y=1:M
  for x=1:N
      
      if mask(y,x) == 0
          continue;
      end

      if y < M  && mask(y+1,x) == 1
          row(cnt) = y+(x-1)*M;
          col(cnt) = y+(x-1)*M;
          val(cnt) = -1;
          cnt = cnt+1;
      
          row(cnt) = y+(x-1)*M;
          col(cnt) = y+(x-1)*M+1;
          val(cnt) = 1;
          cnt = cnt+1;
      else
          row(cnt) = y+(x-1)*M;
          col(cnt) = y+(x-1)*M;
          val(cnt) = 1;
          cnt = cnt+1;
      
          row(cnt) = y+(x-1)*M;
          col(cnt) = y+(x-1)*M-1;
          val(cnt) = -1;
          cnt = cnt+1;
      end
  end
end
row = row(1:cnt-1);
col = col(1:cnt-1);
val = val(1:cnt-1);
Kyf = sparse(row,col,val,M*N,M*N);

nabla = [Kxf; Kyf];


valx = sum(abs(nabla),1) > 0;
nabla_ = nabla(:, valx);
valy = sum(abs(nabla_),2) > 0;
nabla_ = nabla_(valy,:);

% nz = nz(:);
% nz = nz(valx);
% 
% W = spdiags([nz(:); nz(:)], 0, size(nabla_,1), size(nabla_,1));
% 
% A = W*nabla_;