function res = medfilt3(input, dim)

res = zeros(size(input));

parfor y = 1:size(input,1)
    res_temp = zeros(1,size(input,2),size(input,3));
    for x = 1:size(input,2)
        for z = 1:size(input,3)
            minY = max(1, y-dim(1));
            maxY = min(size(input,1), y+dim(1));
            minX = max(1, x-dim(2));
            maxX = min(size(input,2), x+dim(2));
            minZ = max(1, z-dim(3));
            maxZ = min(size(input,3), z+dim(3));
            temp = input(minY:maxY, minX:maxX, minZ:maxZ);
            res_temp(1,x,z) = median(temp(:));
        end
    end
    res(y,:,:) = res_temp;
end
