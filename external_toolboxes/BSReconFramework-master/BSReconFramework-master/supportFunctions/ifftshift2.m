function S = ifftshift2(I)
    S = ifftshift(ifftshift(I,1),2);
end