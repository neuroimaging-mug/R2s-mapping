function pattern = getElipticPattern (Wy, Hz, dataDim)

y0 = floor(dataDim(1)/2);
z0 = floor(dataDim(3)/2);
refRange_y = y0 - floor(Wy/2) : y0 + floor(Wy/2);

pattern = zeros(dataDim(1), dataDim(3));

for i = 1:length(refRange_y)
    zmin = round(z0-sqrt((1-((refRange_y(i)-y0)^2)/((Wy/2)^2))*(Hz/2)^2));
    zmax = round(z0+sqrt((1-((refRange_y(i)-y0)^2)/((Wy/2)^2))*(Hz/2)^2));
    pattern(refRange_y(i), zmin:zmax) = 1;
end

pattern = repmat(pattern, [1,1,dataDim(2)]);
pattern = permute(pattern, [1,3,2]);
pattern = repmat(pattern, [1,1,1,dataDim(4),dataDim(5)]);
