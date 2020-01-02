function K = myIfft(M)

 	si = size(M);
  	a =  sqrt(si(1)) * sqrt(si(2)); 

    K = ifftshift2(ifft2(ifftshift2(M))) * a;
    
    
end
