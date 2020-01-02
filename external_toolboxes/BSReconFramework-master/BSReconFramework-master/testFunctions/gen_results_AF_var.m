imagine(squeeze(pattern_(:,1,:,2,:,:)));



i = 5;
ACL_temp(i)
AF_temp
imagine(flipAngleMap_full, squeeze(flipAngleMap(:,:,:,i,:)))
imagine(squeeze(difference(:,:,:,i,:)))

imagine(squeeze(difference(:,:,:,i,7:9)))
imagine(flipAngleMap_full, squeeze(flipAngleMap(:,:,:,i,4:6)))



%% ROI
y = 30:107;
x =  38:91;
z =   1:22;
ROI = difference(y,x, z, i,:);
ROI = reshape(ROI, size(ROI,1)*size(ROI,2)*size(ROI,3), size(ROI,5));
maxDiff  = max(ROI)
meanDiff = mean(ROI)
stdDiff  = std(ROI)

pattern = pattern_(:,:,:,1,i,10);
numel(pattern)./sum(pattern(:))
imagine(squeeze(pattern(:,1,:)))

figure
boxplot(ROI, 'whisker', 15)
figure
axes('FontSize',14)
boxplot(ROI(:,[1,3,5,8,10]), 'whisker', 15,'labels',{num2str(AF_temp(1)), num2str(AF_temp(3)), num2str(AF_temp(5)), num2str(AF_temp(8)), num2str(AF_temp(10))})
set(findobj(gca,'Type','text'),'FontSize',15)
txt = findobj(gca,'Type','text');
set(txt(1:end),'VerticalAlignment', 'Middle');
title(['ACL-region: n_{ACL} = ', num2str(ACL_temp(i)),'^2'],  'FontSize', 20)
xlabel('AF outside ACL-region', 'FontSize', 16)
ylabel('Absolute difference to reference in %', 'FontSize', 16)
h = get(gca, 'xlabel');
oldpos = get(h, 'Position');
set(h, 'Position', oldpos - [-210, 15, 0])
ax = findobj(gcf,'type','axes');
for j=1:length(ax),
  boxparent = getappdata(ax(j),'boxplothandle');
  listeners = getappdata(boxparent,'boxlisteners');
  for k = 1:length(listeners),
    listeners{k}.Enabled = 'off';
  end
end



temp = flipAngleMap_full.*mask;
figure;
imshow(temp(:,:,16), [80 120])
colormap('jet')
colorbar

figure;imshow(permute(squeeze(temp(64,:,:)),[2 1]), [80 120])
colormap('jet')
colorbar

figure;imshow(squeeze(temp(:,64,:)), [80 120])
colormap('jet')
colorbar