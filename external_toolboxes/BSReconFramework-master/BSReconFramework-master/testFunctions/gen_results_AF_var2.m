imagine(squeeze(pattern_(:,1,:,2,:,:)));




ACL_temp
AF_temp
imagine(flipAngleMap_full, squeeze(flipAngleMap(:,:,:,:,:)))
imagine(squeeze(difference(:,:,:,:,:)))

imagine(squeeze(difference(:,:,:,:,7:9)))
imagine(flipAngleMap_full, squeeze(flipAngleMap(:,:,:,:,10:11)))



%ROI
y = 30:107;
x =  38:91;
z =   1:22;
ROI = difference(y,x, z, :,:);
ROI = reshape(ROI, size(ROI,1)*size(ROI,2)*size(ROI,3), size(ROI,4));
maxDiff  = max(ROI)
meanDiff = mean(ROI)
stdDiff  = std(ROI)

pattern = pattern_(:,:,:,1,5,1);
numel(pattern)./sum(pattern(:))
imagine(squeeze(pattern(:,1,:)))


imagine(flipAngleMap_full, squeeze(flipAngleMap(:,:,:,5)), squeeze(flipAngleMap_lowres(:,:,:,5)))
