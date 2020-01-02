function pattern = getRectPattern (Wy, Hz, dataDim)

refRange = floor(dataDim(1)/2) - floor(Wy/2) + 1 : floor(dataDim(1)/2) + floor(Wy/2);
refRangeSlice = floor(dataDim(3)/2) - floor(Hz/2) + 1 : floor(dataDim(3)/2) + floor(Hz/2);
pattern = zeros(dataDim(1), dataDim(3));
pattern(refRange, refRangeSlice) = 1;
pattern = repmat(pattern, [1,1,dataDim(2)]);
pattern = permute(pattern, [1,3,2]);
pattern = repmat(pattern, [1,1,1,dataDim(4),dataDim(5)]);