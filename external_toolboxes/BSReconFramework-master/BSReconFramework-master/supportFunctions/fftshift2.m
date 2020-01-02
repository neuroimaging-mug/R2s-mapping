function S = fftshift2(I)
    S = fftshift(fftshift(I,1),2);
end