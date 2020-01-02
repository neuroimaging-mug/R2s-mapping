%% ACL 12
temp = flipAngleMap(:,:,:,2).*mask;
figure;
imshow(temp(:,:,16), [75 125])
colormap('jet')
colorbar

figure;imshow(permute(squeeze(temp(64,:,:)),[2 1]), [75 125])
colormap('jet')
colorbar

figure;imshow(squeeze(temp(:,64,:)), [75 125])
colormap('jet')
colorbar

%% reference
temp = flipAngleMap_full.*mask;
figure;
imshow(temp(:,:,16), [75 125])
colormap('jet')
c = colorbar();
set(gca,'FontSize', 13);
zlab = get(c,'ylabel');
set(zlab,'String','rel. B_{1+}-field in % of the nominal flip-angle', 'FontSize', 14); 

figure;imshow(permute(squeeze(temp(64,:,:)),[2 1]), [75 125])
colormap('jet')
colorbar

figure;imshow(squeeze(temp(:,64,:)), [75 125])
colormap('jet')
colorbar
%% difference

y = 30:107;
x =  38:91;
z =   1:22;
temp = difference(:,:,:,5).*mask;
%transversal
figure;
imshow(temp(:,:,16), [0 5])
line([x(1) x(1)],[y(1) y(end)], 'LineWidth',2.5); line([x(1) x(end)],[y(1) y(1)], 'LineWidth',2.5); 
line([x(end) x(1)],[y(end) y(end)], 'LineWidth',2.5); line([x(end) x(end)],[y(end) y(1)], 'LineWidth',2.5);
c = colorbar();
set(gca,'FontSize', 13);
zlab = get(c,'ylabel');
set(zlab,'String','Difference to Reference in %', 'FontSize', 14); 

%coronal
figure;imshow(permute(squeeze(temp(64,:,:)),[2 1]), [0 5])
line([x(1) x(1)],[z(1) z(end)], 'LineWidth',2.5); line([x(1) x(end)],[z(1) z(1)], 'LineWidth',2.5); 
line([x(end) x(1)],[z(end) z(end)], 'LineWidth',2.5); line([x(end) x(end)],[z(end) z(1)], 'LineWidth',2.5); 
colorbar

%sagital
figure;imshow(squeeze(temp(:,64,:)), [0 5])
line([z(1) z(1)],[y(1) y(end)], 'LineWidth',2.5); line([z(1) z(end)],[y(1) y(1)], 'LineWidth',2.5); 
line([z(end) z(1)],[y(end) y(end)], 'LineWidth',2.5); line([z(end) z(end)],[y(end) y(1)], 'LineWidth',2.5);
colorbar


%% Low resolution
temp = flipAngleMap_lowres(:,:,:,7).*mask;
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

%% 
ACL_temp
AF_temp
imagine(flipAngleMap_full, squeeze(flipAngleMap(:,:,:,2)), squeeze(flipAngleMap_lowres(:,:,:,2)));
% imagine(squeeze(difference(:,:,:,:,:)))


%% Errorplot ACL only
difference2 = abs(repmat(flipAngleMap_full, [1,1,1,length(ACL_temp)]) - flipAngleMap_lowres);
y = 30:107;
x =  38:91;
z =   2:22;
ROI = difference(y,x, z, :);
ROI = reshape(ROI, size(ROI,1)*size(ROI,2)*size(ROI,3), size(ROI,4));
ROI2 = difference2(y,x, z, :);
ROI2 = reshape(ROI2, size(ROI2,1)*size(ROI2,2)*size(ROI2,3), size(ROI2,4));
maxDiff  = max(ROI)
meanDiff = mean(ROI)
stdDiff  = std(ROI)
rmsDiff  = rms(ROI)
rmsDiff2  = rms(ROI2)

for k = 1:size(pattern_,5)
    pattern = pattern_(:,:,:,:,k);
    acc(k) = numel(pattern)./sum(pattern(:));
end
acc

%% boxplot(ROI, 'whisker', 15, 'labels', bla)
figure
axes('FontSize',14)
boxplot(ROI, 'whisker', 15,'labels',{[num2str(ACL_temp(1)),'²'], [num2str(ACL_temp(2)),'²'], [num2str(ACL_temp(3)),'²'], [num2str(ACL_temp(4)),'²'], [num2str(ACL_temp(5)),'²'], [num2str(ACL_temp(6)),'²'], [num2str(ACL_temp(7)),'²']})
set(findobj(gca,'Type','text'),'FontSize',15)
txt = findobj(gca,'Type','text');
set(txt(1:end),'VerticalAlignment', 'Middle');
%title('ACL-only',  'FontSize', 20)
grid on
xlabel('ACL-region size: n_{c}', 'FontSize', 16)
ylabel('Deviation to reference in %', 'FontSize', 16)
h = get(gca, 'xlabel');
oldpos = get(h, 'Position');
set(h, 'Position', oldpos - [0, 15, 0])
ax = findobj(gcf,'type','axes');
for j=1:length(ax),
  boxparent = getappdata(ax(j),'boxplothandle');
  listeners = getappdata(boxparent,'boxlisteners');
  for k = 1:length(listeners),
    listeners{k}.Enabled = 'off';
  end
end