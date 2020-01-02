
n=10;
m=150;
k=12;
nch = 5;


dx = 0.1;
dy = 2;
dz = 1.3;


%Vector operators
x = rand(n,m,k,nch);
y = rand(n,m,k,3,nch);

gx = fgrad_3(x,dx,dy,dz);
divy = bdiv_3(y,dx,dy,dz);

s1 = gx(:)'*y(:);
s2 = -x(:)'*divy(:);

display('############## Vector operators ##############')
display(['S1:      ',num2str(s1)]);
display(['S2:      ',num2str(s2)]);
display(['Dif:     ',num2str(abs(s1-s2))]);
display(['Rel-Dif: ',num2str(abs(s1-s2)/(n*m*k*nch))]);




%Matrix operators
x = rand(n,m,k,3,nch);
y = rand(n,m,k,6,nch);



gx = sym_bgrad_3(x,dx,dy,dz);
divy = fdiv_3(y,dx,dy,dz);



tmp = gx(:,:,:,1:3,:).*y(:,:,:,1:3,:) + 2*gx(:,:,:,4:6,:).*y(:,:,:,4:6,:);

s1 = sum(tmp(:));
s2 = -x(:)'*divy(:);

display('############## Matrix operators ##############')
display(['S1:      ',num2str(s1)]);
display(['S2:      ',num2str(s2)]);
display(['Dif:     ',num2str(abs(s1-s2))]);
display(['Rel-Dif: ',num2str(abs(s1-s2)/(n*m*k*nch))]);


%% MR OP

n=10;
m=150;
k=12;
nch = 5;


x = rand(n,m,k);
y = rand(n,m,k,nch);

mri_obj_test.b1 = ones(n,m,k,nch);
mri_obj_test.mask = ones(n,m,k);

KHy = backward_mri3d(y, mri_obj_test.b1, mri_obj_test.mask);

Kx = forward_mri3d(x, mri_obj_test.b1, mri_obj_test.mask);

s1 = dot(y(:),Kx(:));
s2 = dot(x(:),KHy(:));


display('############## Matrix operators ##############')
display(['S1:      ',num2str(s1)]);
display(['S2:      ',num2str(s2)]);
display(['Dif:     ',num2str(abs(s1-s2))]);
display(['Rel-Dif: ',num2str(abs(s1-s2)/(n*m*k))]);
%%

nn = 100;
mm = 101;
kk = 101;

x = rand(nn,mm,kk);
y = rand(nn,mm,kk);

%Ax = fftn(x)./sqrt(nn*mm*kk);
%AHy = ifftn(x).*sqrt(nn*mm*kk);

Ax = fft3c(x);
AHy = ifft3c(x);

s1 = dot(y(:),Ax(:));
s2 = dot(x(:),AHy(:));


display('############## Matrix operators ##############')
display(['S1:      ',num2str(s1)]);
display(['S2:      ',num2str(s2)]);
display(['Dif:     ',num2str(abs(s1-s2))]);
display(['Rel-Dif: ',num2str(abs(s1-s2)/(nn*mm*kk))]);

%%
[n,m,l] = size(testdata);
[x,y,z] = meshgrid(1:m,1:n,1:l);
chop = (-1).^(x+y+z+1);
clear x y z
testdata2 = fftshift(chop.*testdata);


image1 = ifft3c(testdata);
image2 = ifftn(testdata2).*sqrt(n*m*l);
image3 = image_shift3d(fftn(testdata2)./sqrt(n*m*l));

imagine(image1,image2,image3)


a = (1:10);

b1 = ifft(fft(a));
b2 = fft(ifft(a));

b3 = ifft(a).*sqrt(length(a));
b4 = fft(a)./sqrt(length(a));
b5 = flip(b4);
b6 = real(b4)-1i.*imag(b4);

figure;
subplot(1,2,1); hold on; plot(real(b3),'-xr'); plot(real(b6),'-sr');
subplot(1,2,2); hold on; plot(imag(b3),'-xb'); plot(imag(b4),'-ob');plot(imag(b5),'-db');plot(imag(b6),'-sb')

%%


a = meshgrid(1:10,1:10);
a = ifft2(a);

b1 = ifft2(fft2(a));
b2 = fft2(ifft2(a));

idx = b1~=b2;
sum(b1(:)==b2(:))./numel(b1)

b3 = ifft2(a).*sqrt(length(a(:)));
b4 = fft2(a)./sqrt(length(a(:)));


b5 = flip(flip(circshift(circshift(b4,-1,1),-1,2),1),2);

sum(b3(:)==b5(:))./numel(b3)
idx = b3~=b5;



