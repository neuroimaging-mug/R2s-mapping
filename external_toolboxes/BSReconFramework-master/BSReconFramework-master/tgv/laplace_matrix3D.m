
function A = laplace_matrix3D(n,m,z,u,bcond)


A = sparse(n*m*z,n*m*z);
if ~isempty(u)
    u = u(:);
end


%Inner diagonal cubes

	for j=1:n*m*z
		A(j,j) = 6; %Set diagonals
	end
    
	for j=1:n*m*z -1
		A(j+1,j) = -1;
		A(j,j+1) = -1;
	end

%Off diagonal cubes (2nd dimension)
	
    for j=1:n*(m-1)
        for k=1:z
            A((k-1)*m*n+n+j,(k-1)*m*n+j) = -1;
            A((k-1)*m*n+j,(k-1)*m*n+n+j) = -1;
        end
    end
    
%Off diagonal cubes (3rd dimension)
	
	for j=1:n*m*(z-1)
		A(n*m+j,j) = -1;
		A(j,n*m+j) = -1;
	end

%Periodic boundary extension
if strcmp(bcond,'per')
	
	for j=n:n:n*m*z -1
		A(j+1,j) = 0;
		A(j,j+1) = 0;
    end
    for j=n*m:n*m:n*m*z -1
		A(j+n,j) = 0;
		A(j,j+n) = 0;
	end
	
	
	for j=1:n %Top right and lower left block (periodic row)
        for k=1:z
            A((k-1)*m*n+j,(k-1)*m*n+n*(m-1)+j) = -1;
            A((k-1)*m*n+n*(m-1)+j,(k-1)*m*n+j) = -1;
        end
	end
	
	for j=1:n:n*m %Peroidic column
        for k=1:z
            A((k-1)*m*n+j,(k-1)*m*n+j+(n-1)) = -1;
            A((k-1)*m*n+j+(n-1),(k-1)*m*n+j) = -1;
        end
    end
    
    for j=1:n*m %Top right and lower left block (periodic slice)
		A(j,n*m*(z-1)+j) = -1;
		A(n*m*(z-1)+j,j) = -1;
	end
	

%Symmetric boundary extension
elseif strcmp(bcond,'sym')


	for j=1:n:n*m - 1 %Set boundary of inner diagonal cubes
		A(j,j) = 3;
		A(j + n-1, j+n-1) = 3;
	end

	for j=1:n %Reduce first and last diagonal cube
		A(j,j) = A(j,j) - 1;
		A(n*(m - 1) + j,n*(m - 1) + j) = A(n*(m - 1) + j,n*(m - 1) + j) - 1;
	end



	for j=n:n:n*m -1
		A(j+1,j) = 0;
		A(j,j+1) = 0;
	end
	
end	


%Add diagonal entries
if ~isempty(u)
        for j=1:n*m
            A(j,j) = A(j,j) + u(j);
        end
end
