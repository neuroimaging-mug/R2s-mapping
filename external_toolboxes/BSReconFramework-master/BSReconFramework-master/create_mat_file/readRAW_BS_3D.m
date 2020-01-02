function [data, hdr, PhaseEncDir] = readRAW_BS_3D(inFile)

image_obj = mapVBVD(inFile,'removeOS');

if length(image_obj) == 2
    image_obj = image_obj{2};
end

dataSize = length(image_obj.image.sqzSize);

hdr = image_obj.hdr;
if ~isfield(hdr.MeasYaps.sSliceArray.asSlice{1}, 'dInPlaneRot')
    PhaseEncDir = 'LIN';
%     
%     %MS: try to change PhaseEncDir to col
%       PhaseEncDir = 'COL';
else
    PhaseEncDir = 'COL';
end

MRAcquisitionType = hdr.Dicom.tMRAcquisitionType;

switch dataSize
    case 4
        data = image_obj.image{:,:,:,:};
        
        switch PhaseEncDir
            case 'LIN'
                data = permute(data, [3,1,4,2]);
            case 'COL'
                data = permute(data, [1,3,4,2]);
        end
        
        data = flip(data, 1);
        data = reshape(data, [size(data,1), size(data,2), 1, size(data,3), size(data,4)]);
        
    case 5
        data = image_obj.image{:,:,:,:,:};
        
        switch PhaseEncDir
            case 'LIN'
                data = permute(data, [3,1,4,5,2]);
            case 'COL'
                data = permute(data, [1,3,4,5,2]);
        end
        
        data = flip(data, 1);
        
        if strcmp(MRAcquisitionType, '3D')
            data = flip(data, 3);
        end
        
%         if strcmp(MRAcquisitionType, '3D')
%             data = flip(data, 3);
%             temp = ifftshift(ifft(ifftshift(data,3), [], 3),3);
%             sliceDiff = (image_obj.image.NPar - numSlice)/2;
%             temp = temp(:,:,sliceDiff+1:end-sliceDiff,:,:);
%             data = fftshift(fft(fftshift(temp,3), [], 3),3);
%         end
        
    otherwise
        error('Dimension not supported');

end
end





