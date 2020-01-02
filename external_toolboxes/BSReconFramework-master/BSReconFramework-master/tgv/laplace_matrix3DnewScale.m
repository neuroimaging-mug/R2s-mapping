
function A = laplace_matrix3DnewScale(n,m,z,u,bcond, dz)

%A = sparse(n*m*z,n*m*z);
if ~isempty(u)
    u = u(:);
end
i = []; j = []; s = []; 
nmz = n*m*z;
dz = 1/(dz^2);

%Periodic boundary extension
if strcmp(bcond,'per')
    
%Inner diagonal cubes
    %Set diagonals
    i1 = 1:nmz;
    j1 = 1:nmz;
    s1 = 6*ones(1, length(i1));
    i = [i, i1]; j = [j, j1]; s = [s, s1];
    
    i1 = repmat(repmat((1:n-1)+1, [m, 1]) + repmat(((1:m)'-1)*n, [1, n-1]), [1, 1, z]) + repmat(reshape(((1:z)-1) * n*m,[1,1,z]), [m, n-1, 1]);
    j1 = repmat(repmat((1:n-1), [m, 1]) + repmat(((1:m)'-1)*n, [1, n-1]), [1, 1, z]) + repmat(reshape(((1:z)-1) * n*m,[1,1,z]), [m, n-1, 1]);
    i1 = i1(:)';
    j1 = j1(:)';
    s1 = -ones(1, length(i1)+length(j1));
    i = [i, i1, j1]; j = [j, j1, i1]; s = [s, s1];

%Off diagonal cubes (2nd dimension)
    lauf1=1:n*(m-1);
    lauf2=(1:z);
    i1 = repmat((lauf2-1)*m*n + n, [length(lauf1), 1])+repmat(lauf1', [1, length(lauf2)]);
    i1 = i1(:)'; 
    j1 = repmat((lauf2-1)*m*n, [length(lauf1), 1]) + repmat(lauf1', [1, length(lauf2)]);
    j1 = j1(:)';
    s1 = -ones(1, length(i1)+length(j1));
    i = [i, i1, j1]; j = [j, j1, i1]; s = [s, s1];
    
%Off diagonal cubes (3rd dimension)
    i1 = (1:n*m*(z-1)) + n*m;
    j1 = 1:n*m*(z-1);
    s1 = -ones(1, length(i1)+length(j1));
    i = [i, i1, j1]; j = [j, j1, i1]; s = [s, s1];
	
    %Top right and lower left block (periodic row)
    lauf1=1:n;
    lauf2=1:z;
    i1 = repmat((lauf2-1)*m*n, [length(lauf1),1])+repmat(lauf1', [1, length(lauf2)]);
    j1 = repmat((lauf2-1)*m*n+n*(m-1), [length(lauf1),1])+repmat(lauf1', [1, length(lauf2)]);
    i1 = i1(:)';
    j1 = j1(:)';
    s1 = -ones(1, length(i1)+length(j1));
    i = [i, i1, j1]; j = [j, j1, i1]; s = [s, s1];
    
    %Peroidic column

    lauf1=1:n:n*m;
    lauf2=1:z;
    i1 = repmat((lauf2-1)*m*n, [length(lauf1),1])+repmat(lauf1', [1, length(lauf2)]);
    j1 = repmat((lauf2-1)*m*n+(n-1), [length(lauf1),1])+repmat(lauf1', [1, length(lauf2)]);
    i1 = i1(:)';
    j1 = j1(:)';
    s1 = -ones(1, length(i1)+length(j1));
    i = [i, i1, j1]; j = [j, j1, i1]; s = [s, s1];
    
    %Top right and lower left block (periodic slice)
    i1 = 1:n*m;
    j1 = n*m*(z-1)+(1:n*m);
    s1 = -ones(1, length(i1)+length(j1));
    i = [i, i1, j1]; j = [j, j1, i1]; s = [s, s1];
	

%Symmetric boundary extension
elseif strcmp(bcond,'sym')

    %Inner diagonal cubes
    %Set diagonals
    i1 = 1:nmz;
    j1 = 1:nmz;
%     s1 = 6*ones(1, length(i1));
    s1 = (4+2*dz)*ones(1, length(i1));
    i = [i, i1]; j = [j, j1]; s = [s, s1];
    
    nm = n*m;
    % ---------------------- Outer corners val 1st and last subdiagonal element 3 ----------------------
    % 1st dim
    i1 = [1,n, n*m, n*m-n+1, nmz, nmz-n*m+1, nmz-n+1, nmz-n*m+n];
    j1 = [1+1,n-1, n*m-1, n*m-n+1+1, nmz-1, nmz-n*m+1+1, nmz-n+1+1, nmz-n*m+n-1];
    s1 = -2*ones(1, length(i1));
    i = [i, i1]; j = [j, j1]; s = [s, s1];
    % 2nd dim
    i1 = [1,n, n*m, n*m-n+1, nmz, nmz-n*m+1, nmz-n+1, nmz-n*m+n];
    j1 = [1+n,n+n, n*m-n, n*m-n+1-n, nmz-n, nmz-n*m+1+n, nmz-n+1-n, nmz-n*m+n+n];
    s1 = -2*ones(1, length(i1));
    i = [i, i1]; j = [j, j1]; s = [s, s1];
    % 3rd dim
    i1 = [1,n, n*m, n*m-n+1, nmz, nmz-n*m+1, nmz-n+1, nmz-n*m+n];
    j1 = [1+nm,n+nm, n*m+nm, n*m-n+1+nm, nmz-nm, nmz-n*m+1-nm, nmz-n+1-nm, nmz-n*m+n-nm];
    s1 = -2*dz*ones(1, length(i1));
    i = [i, i1]; j = [j, j1]; s = [s, s1];
    
    
    % ------------------ Outer Elements elements 1st and last subdiagonal element val 4 ------------------------
    % 1st dim
    i1 = [2:n-1, n*m-1:-1:n*m-n+2, nmz-1:-1:nmz-n+2, nmz-n*m+2:nmz-n*m+n-1];
    j1 = i1 + 1;
    s1 = -1*ones(1, length(i1));
    i = [i, i1]; j = [j, j1]; s = [s, s1];
    
    i1 = [2:n-1, n*m-1:-1:n*m-n+2, nmz-1:-1:nmz-n+2, nmz-n*m+2:nmz-n*m+n-1];
    j1 = i1 - 1;
    s1 = -1*ones(1, length(i1));
    i = [i, i1]; j = [j, j1]; s = [s, s1];
    % 2nd dim
    i1 = [2:n-1, n*m-1:-1:n*m-n+2, nmz-1:-1:nmz-n+2, nmz-n*m+2:nmz-n*m+n-1];
    j1 = [(2:n-1)+n, (n*m-1:-1:n*m-n+2)-n, (nmz-1:-1:nmz-n+2)-n, (nmz-n*m+2:nmz-n*m+n-1)+n];
    s1 = -2*ones(1, length(i1));
    i = [i, i1]; j = [j, j1]; s = [s, s1];
    % 3rd dim
    i1 = [2:n-1, n*m-1:-1:n*m-n+2, nmz-1:-1:nmz-n+2, nmz-n*m+2:nmz-n*m+n-1];
    j1 = [(2:n-1)+nm, (n*m-1:-1:n*m-n+2)+nm, (nmz-1:-1:nmz-n+2)-nm, (nmz-n*m+2:nmz-n*m+n-1)-nm];
    s1 = -2*dz*ones(1, length(i1));
    i = [i, i1]; j = [j, j1]; s = [s, s1];
    
    
    
    
     % ------------------ Inner Elements elements 1st and last subdiagonal element val 5 ------------------------
     % 1st dim
    lauf1 = [2:n-1, nmz-n*m+2:nmz-n*m+n-1]; lauf2 = (2:m-1)-1;
    i1 = repmat(lauf1, [length(lauf2),1]) + n*repmat(lauf2', [1, length(lauf1)]);
    i1 = i1(:)';
    j1 = i1+1;
    s1 = -1*ones(1, length(i1));
    i = [i, i1]; j = [j, j1]; s = [s, s1];
    
    lauf1 = [2:n-1, nmz-n*m+2:nmz-n*m+n-1]; lauf2 = (2:m-1)-1;
    i1 = repmat(lauf1, [length(lauf2),1]) + n*repmat(lauf2', [1, length(lauf1)]);
    i1 = i1(:)';
    j1 = i1-1;
    s1 = -1*ones(1, length(i1));
    i = [i, i1]; j = [j, j1]; s = [s, s1];
    
    %2nd dim
    lauf1 = [2:n-1, nmz-n*m+2:nmz-n*m+n-1]; lauf2 = (2:m-1)-1;
    i1 = repmat(lauf1, [length(lauf2),1]) + n*repmat(lauf2', [1, length(lauf1)]);
    i1 = i1(:)';
    j1 = i1+n;
    s1 = -1*ones(1, length(i1));
    i = [i, i1]; j = [j, j1]; s = [s, s1];
    
    lauf1 = [2:n-1, nmz-n*m+2:nmz-n*m+n-1]; lauf2 = (2:m-1)-1;
    i1 = repmat(lauf1, [length(lauf2),1]) + n*repmat(lauf2', [1, length(lauf1)]);
    i1 = i1(:)';
    j1 = i1-n;
    s1 = -1*ones(1, length(i1));
    i = [i, i1]; j = [j, j1]; s = [s, s1];
    
    %3rd dim
    lauf1 = [2:n-1, nmz-n*m+2:nmz-n*m+n-1]; lauf2 = (2:m-1)-1;
    i1 = repmat(lauf1, [length(lauf2),1]) + n*repmat(lauf2', [1, length(lauf1)]);
    i1 = i1(:)';
    lauf3 = [(2:n-1)+nm, (nmz-n*m+2:nmz-n*m+n-1)-nm]; lauf4 = (2:m-1)-1;
    j1 = repmat(lauf3, [length(lauf4),1]) + n*repmat(lauf4', [1, length(lauf3)]);
    j1 = j1(:)';
    s1 = -2*dz*ones(1, length(i1));
    i = [i, i1]; j = [j, j1]; s = [s, s1];
    
    
    % ------------------ Inner Elements corner elements 1st and last subdiagonal element val 4 ------------------------
    % 1st dim 
    lauf1 = [1, n, nmz-n*m+1, nmz-n*m+n]; lauf2 = (2:m-1)-1;
    i1 = repmat(lauf1, [length(lauf2),1]) + n*repmat(lauf2', [1, length(lauf1)]);
    i1 = i1(:)';
    lauf3 = [1+1, n-1, nmz-n*m+1+1, nmz-n*m+n-1]; lauf4 = (2:m-1)-1;
    j1 = repmat(lauf3, [length(lauf4),1]) + n*repmat(lauf4', [1, length(lauf3)]);
    j1 = j1(:)';
    s1 = -2*ones(1, length(i1));
    i = [i, i1]; j = [j, j1]; s = [s, s1];
    
    % 2nd dim 
    lauf1 = [1, n, nmz-n*m+1, nmz-n*m+n]; lauf2 = (2:m-1)-1;
    i1 = repmat(lauf1, [length(lauf2),1]) + n*repmat(lauf2', [1, length(lauf1)]);
    i1 = i1(:)';
    j1 = i1+n;
    s1 = -1*ones(1, length(i1));
    i = [i, i1]; j = [j, j1]; s = [s, s1];
    
    lauf1 = [1, n, nmz-n*m+1, nmz-n*m+n]; lauf2 = (2:m-1)-1;
    i1 = repmat(lauf1, [length(lauf2),1]) + n*repmat(lauf2', [1, length(lauf1)]);
    i1 = i1(:)';
    j1 = i1-n;
    s1 = -1*ones(1, length(i1));
    i = [i, i1]; j = [j, j1]; s = [s, s1];
    
    % 3rd dim 
    lauf1 = [1, n, nmz-n*m+1, nmz-n*m+n]; lauf2 = (2:m-1)-1;
    i1 = repmat(lauf1, [length(lauf2),1]) + n*repmat(lauf2', [1, length(lauf1)]);
    i1 = i1(:)';
    lauf3 = [1+nm, n+nm, nmz-n*m+1-nm, nmz-n*m+n-nm]; lauf4 = (2:m-1)-1;
    j1 = repmat(lauf3, [length(lauf4),1]) + n*repmat(lauf4', [1, length(lauf3)]);
    j1 = j1(:)';
    s1 = -2*dz*ones(1, length(i1));
    i = [i, i1]; j = [j, j1]; s = [s, s1];
    

    
    % ---------------------- Outer corners inner subdiagonal element val 4 ----------------------
    % 1st dim
    lauf1 = [1,n, n*m, n*m-n+1]; lauf2 = (2:z-1)-1;
    i1 = repmat(lauf1, [length(lauf2),1]) + n*m*repmat(lauf2', [1, length(lauf1)]);
    i1 = i1(:)';
    lauf3 = [1+1,n-1, n*m-1, n*m-n+1+1]; lauf4 = (2:z-1)-1;
    j1 = repmat(lauf3, [length(lauf4),1]) + n*m*repmat(lauf4', [1, length(lauf3)]);
    j1 = j1(:)';
    s1 = -2*ones(1, length(i1));
    i = [i, i1]; j = [j, j1]; s = [s, s1];
    % 2nd dim
    lauf1 = [1,n, n*m, n*m-n+1]; lauf2 = (2:z-1)-1;
    i1 = repmat(lauf1, [length(lauf2),1]) + n*m*repmat(lauf2', [1, length(lauf1)]);
    i1 = i1(:)';
    lauf3 = [1+n,n+n, n*m-n, n*m-n+1-n]; lauf4 = (2:z-1)-1;
    j1 = repmat(lauf3, [length(lauf4),1]) + n*m*repmat(lauf4', [1, length(lauf3)]);
    j1 = j1(:)';
    s1 = -2*ones(1, length(i1));
    i = [i, i1]; j = [j, j1]; s = [s, s1];
    % 3rd dim
    lauf1 = [1,n, n*m, n*m-n+1]; lauf2 = (2:z-1)-1;
    i1 = repmat(lauf1, [length(lauf2),1]) + n*m*repmat(lauf2', [1, length(lauf1)]);
    i1 = i1(:)';
    j1 = i1+nm;
    s1 = -1*dz*ones(1, length(i1));
    i = [i, i1]; j = [j, j1]; s = [s, s1];
    
    lauf1 = [1,n, n*m, n*m-n+1]; lauf2 = (2:z-1)-1;
    i1 = repmat(lauf1, [length(lauf2),1]) + n*m*repmat(lauf2', [1, length(lauf1)]);
    i1 = i1(:)';
    j1 = i1-nm;
    s1 = -1*dz*ones(1, length(i1));
    i = [i, i1]; j = [j, j1]; s = [s, s1];
   
    
    % ------------------ Outer Elements  inner subdiagonal elements val 5 ------------------------
    % 1st dim
    lauf1 = [2:n-1, n*m-1:-1:n*m-n+2]; lauf2 = (2:z-1)-1;
    i1 = repmat(lauf1, [length(lauf2),1]) + n*m*repmat(lauf2', [1, length(lauf1)]);
    i1 = i1(:)';
    j1 = i1+1;
    s1 = -1*ones(1, length(i1));
    i = [i, i1]; j = [j, j1]; s = [s, s1];
    
    lauf1 = [2:n-1, n*m-1:-1:n*m-n+2]; lauf2 = (2:z-1)-1;
    i1 = repmat(lauf1, [length(lauf2),1]) + n*m*repmat(lauf2', [1, length(lauf1)]);
    i1 = i1(:)';
    j1 = i1-1;
    s1 = -1*ones(1, length(i1));
    i = [i, i1]; j = [j, j1]; s = [s, s1];
    
    % 2nd dim
    lauf1 = [2:n-1, n*m-1:-1:n*m-n+2]; lauf2 = (2:z-1)-1;
    i1 = repmat(lauf1, [length(lauf2),1]) + n*m*repmat(lauf2', [1, length(lauf1)]);
    i1 = i1(:)';
    lauf3 = [(2:n-1)+n, (n*m-1:-1:n*m-n+2)-n]; lauf4 = (2:z-1)-1;
    j1 = repmat(lauf3, [length(lauf4),1]) + n*m*repmat(lauf4', [1, length(lauf3)]);
    j1 = j1(:)';
    s1 = -2*ones(1, length(i1));
    i = [i, i1]; j = [j, j1]; s = [s, s1];
    
    % 3rd dim
    lauf1 = [2:n-1, n*m-1:-1:n*m-n+2]; lauf2 = (2:z-1)-1;
    i1 = repmat(lauf1, [length(lauf2),1]) + n*m*repmat(lauf2', [1, length(lauf1)]);
    i1 = i1(:)';
    j1 = i1+nm;
    s1 = -1*dz*ones(1, length(i1));
    i = [i, i1]; j = [j, j1]; s = [s, s1];
    
    lauf1 = [2:n-1, n*m-1:-1:n*m-n+2]; lauf2 = (2:z-1)-1;
    i1 = repmat(lauf1, [length(lauf2),1]) + n*m*repmat(lauf2', [1, length(lauf1)]);
    i1 = i1(:)';
    j1 = i1-nm;
    s1 = -1*dz*ones(1, length(i1));
    i = [i, i1]; j = [j, j1]; s = [s, s1];
    
    
    % ---------------------------- Inner elements inner subdiagonal elements val 6 ----------------------------------
    % 1st dim
    lauf1 = 2:n-1; lauf2 = (2:m-1)-1; lauf3 = (2:z-1)-1;
    i1 = (repmat(lauf1, [length(lauf2),1,length(lauf3)]) + n*repmat(lauf2', [1, length(lauf1),length(lauf3)])) + ...
        n*m*repmat(reshape(lauf3,[1,1,length(lauf3)]), [length(lauf2), length(lauf1), 1]);
    i1 = i1(:)';
    j1 = i1+1;
    s1 = -1*ones(1, length(i1));
    i = [i, i1]; j = [j, j1]; s = [s, s1];
    
    lauf1 = 2:n-1; lauf2 = (2:m-1)-1; lauf3 = (2:z-1)-1;
    i1 = (repmat(lauf1, [length(lauf2),1,length(lauf3)]) + n*repmat(lauf2', [1, length(lauf1),length(lauf3)])) + ...
        n*m*repmat(reshape(lauf3,[1,1,length(lauf3)]), [length(lauf2), length(lauf1), 1]);
    i1 = i1(:)';
    j1 = i1-1;
    s1 = -1*ones(1, length(i1));
    i = [i, i1]; j = [j, j1]; s = [s, s1];
    
    % 2nd dim
    lauf1 = 2:n-1; lauf2 = (2:m-1)-1; lauf3 = (2:z-1)-1;
    i1 = (repmat(lauf1, [length(lauf2),1,length(lauf3)]) + n*repmat(lauf2', [1, length(lauf1),length(lauf3)])) + ...
        n*m*repmat(reshape(lauf3,[1,1,length(lauf3)]), [length(lauf2), length(lauf1), 1]);
    i1 = i1(:)';
    j1 = i1+n;
    s1 = -1*ones(1, length(i1));
    i = [i, i1]; j = [j, j1]; s = [s, s1];
    
    lauf1 = 2:n-1; lauf2 = (2:m-1)-1; lauf3 = (2:z-1)-1;
    i1 = (repmat(lauf1, [length(lauf2),1,length(lauf3)]) + n*repmat(lauf2', [1, length(lauf1),length(lauf3)])) + ...
        n*m*repmat(reshape(lauf3,[1,1,length(lauf3)]), [length(lauf2), length(lauf1), 1]);
    i1 = i1(:)';
    j1 = i1-n;
    s1 = -1*ones(1, length(i1));
    i = [i, i1]; j = [j, j1]; s = [s, s1];
    
    % 3rd dim
    lauf1 = 2:n-1; lauf2 = (2:m-1)-1; lauf3 = (2:z-1)-1;
    i1 = (repmat(lauf1, [length(lauf2),1,length(lauf3)]) + n*repmat(lauf2', [1, length(lauf1),length(lauf3)])) + ...
        n*m*repmat(reshape(lauf3,[1,1,length(lauf3)]), [length(lauf2), length(lauf1), 1]);
    i1 = i1(:)';
    j1 = i1+nm;
    s1 = -1*dz*ones(1, length(i1));
    i = [i, i1]; j = [j, j1]; s = [s, s1];
    
    lauf1 = 2:n-1; lauf2 = (2:m-1)-1; lauf3 = (2:z-1)-1;
    i1 = (repmat(lauf1, [length(lauf2),1,length(lauf3)]) + n*repmat(lauf2', [1, length(lauf1),length(lauf3)])) + ...
        n*m*repmat(reshape(lauf3,[1,1,length(lauf3)]), [length(lauf2), length(lauf1), 1]);
    i1 = i1(:)';
    j1 = i1-nm;
    s1 = -1*dz*ones(1, length(i1));
    i = [i, i1]; j = [j, j1]; s = [s, s1];
    
    % -------------------------- val 5 ------------------------------------
    % 1st dim
    lauf1 = [1, n]; lauf2 = (2:m-1)-1; lauf3 = (2:z-1)-1;
    i1 = (repmat(lauf1, [length(lauf2),1,length(lauf3)]) + n*repmat(lauf2', [1, length(lauf1),length(lauf3)])) + ...
        n*m*repmat(reshape(lauf3,[1,1,length(lauf3)]), [length(lauf2), length(lauf1), 1]);
    i1 = i1(:)';
    lauf1 = [1+1, n-1]; lauf2 = (2:m-1)-1; lauf3 = (2:z-1)-1;
    j1 = (repmat(lauf1, [length(lauf2),1,length(lauf3)]) + n*repmat(lauf2', [1, length(lauf1),length(lauf3)])) + ...
        n*m*repmat(reshape(lauf3,[1,1,length(lauf3)]), [length(lauf2), length(lauf1), 1]);
    j1 = j1(:)';
    s1 = -2*ones(1, length(i1));
    i = [i, i1]; j = [j, j1]; s = [s, s1];
    
    % 2nd dim
    lauf1 = [1, n]; lauf2 = (2:m-1)-1; lauf3 = (2:z-1)-1;
    i1 = (repmat(lauf1, [length(lauf2),1,length(lauf3)]) + n*repmat(lauf2', [1, length(lauf1),length(lauf3)])) + ...
        n*m*repmat(reshape(lauf3,[1,1,length(lauf3)]), [length(lauf2), length(lauf1), 1]);
    i1 = i1(:)';
    j1 = i1 + n;
    s1 = -1*ones(1, length(i1));
    i = [i, i1]; j = [j, j1]; s = [s, s1];
    
    lauf1 = [1, n]; lauf2 = (2:m-1)-1; lauf3 = (2:z-1)-1;
    i1 = (repmat(lauf1, [length(lauf2),1,length(lauf3)]) + n*repmat(lauf2', [1, length(lauf1),length(lauf3)])) + ...
        n*m*repmat(reshape(lauf3,[1,1,length(lauf3)]), [length(lauf2), length(lauf1), 1]);
    i1 = i1(:)';
    j1 = i1 - n;
    s1 = -1*ones(1, length(i1));
    i = [i, i1]; j = [j, j1]; s = [s, s1];
    
    % 3rd dim
    lauf1 = [1, n]; lauf2 = (2:m-1)-1; lauf3 = (2:z-1)-1;
    i1 = (repmat(lauf1, [length(lauf2),1,length(lauf3)]) + n*repmat(lauf2', [1, length(lauf1),length(lauf3)])) + ...
        n*m*repmat(reshape(lauf3,[1,1,length(lauf3)]), [length(lauf2), length(lauf1), 1]);
    i1 = i1(:)';
    j1 = i1+nm;
    s1 = -1*dz*ones(1, length(i1));
    i = [i, i1]; j = [j, j1]; s = [s, s1];
    
    lauf1 = [1, n]; lauf2 = (2:m-1)-1; lauf3 = (2:z-1)-1;
    i1 = (repmat(lauf1, [length(lauf2),1,length(lauf3)]) + n*repmat(lauf2', [1, length(lauf1),length(lauf3)])) + ...
        n*m*repmat(reshape(lauf3,[1,1,length(lauf3)]), [length(lauf2), length(lauf1), 1]);
    i1 = i1(:)';
    j1 = i1-nm;
    s1 = -1*dz*ones(1, length(i1));
    i = [i, i1]; j = [j, j1]; s = [s, s1];
end	

A = sparse(i,j,s,nmz, nmz);


%Add diagonal entries
if ~isempty(u)
    i = 1:nmz;
    j = 1:nmz;
    U = sparse(i,j,u,nmz,nmz);
    A = A + U;
end
