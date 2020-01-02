function img = sumSqrImg(img,dim)
% usage: img = sumSqrImg(img,dim)

if nargin ==1
    dim = length(size(img)); % use the last dimensions
end

img = sqrt(sum(img .* conj(img),dim));
end