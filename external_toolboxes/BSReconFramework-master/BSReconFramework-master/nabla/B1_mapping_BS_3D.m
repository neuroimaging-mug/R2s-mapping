function [flipAngleMap, B1Map] = B1_mapping_BS_3D(PathName, FileName1, pulse, deltaf, alphaBS, dims)

dimY = dims(1);
dimX = dims(2);
slices = dims(3);

flipAngleMap = zeros(dimY, dimX, slices);
B1Map = zeros(dimY, dimX, slices);
start_idx = str2double(FileName1(19:20));

parfor i = 1:slices
    idx = start_idx+i-1;
    if idx < 10
        FileName = [FileName1(1:19), num2str(idx), FileName1(21:end)];
    else
        FileName = [FileName1(1:18), num2str(idx), FileName1(21:end)];
    end    
    [flipAngleMap(:,:,i), B1Map(:,:,i)] = B1_mapping_BS(PathName, FileName, pulse, deltaf, alphaBS);
end
end