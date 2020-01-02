function K = myFft(M)

 	si = size(M);
	a = 1 / (sqrt(si(1)) * sqrt(si(2)));      
            
    K = fftshift2(fft2(fftshift2(M))) * a;    
        

end
