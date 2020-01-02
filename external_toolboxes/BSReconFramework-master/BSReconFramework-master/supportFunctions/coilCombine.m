function img = coilCombine(multicoilImgs, C, weightWithProfiles)
% img = coilCombine(multicoilImgs, C, weightWithProfiles)

if nargin < 3
    weightWithProfiles = false;
end

if ~isa(multicoilImgs,'double');
    warning('myRhoEstMultiecho:datatype','converting datatype to double...!');
    multicoilImgs = double(multicoilImgs);
end

[dimY, dimX, NE, NC] = size(multicoilImgs);


img = zeros([dimY,dimX,NE]);

sumC = 1./sum(abs(C).^2,4);

for i = 1 : NE
     img(:,:,i) = sum(ftimes((multicoilImgs(:,:,i,:)) .* conj(C) , sumC),4);
end

if weightWithProfiles
    if NC > 1
        sumSqrCoils = sumSqrImg(C);
    else
        sumSqrCoils = abs(C);
    end
    for i = 1 : NE
        img(:,:,i) = img(:,:,i) .* sumSqrCoils;
    end
end
end