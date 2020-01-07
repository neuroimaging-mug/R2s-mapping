function [ mask ] = thresholdImage(img,th )
%THRESHOLDIMAGE Returns a binary mask for the input img with the selected
%threshold vector th = [th_low, th_up]; 

mask_low = img;  
mask_low(img < th(1)) = 0; 
mask_low(img >= th(1)) = 1; 

mask_high = ones(size(img));  
mask_high(img >= th(2)) = 0; 
mask = mask_low.*mask_high; 

end

