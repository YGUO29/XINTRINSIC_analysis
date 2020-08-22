Max = max(abs(DataMat_mean(:)));
figure, imagesc(DataMat_mean, [-Max, Max]), axis image, colormap(jet), colorbar
h = images.roi.Circle(gca,'Center',[floor(para.width/2) floor(para.height/2)],'Radius',floor(para.height/2)); 
title('Press Enter after the position is adjusted')
pause
mask = createMask(h);
%%
figure,
DataMat_test1 = squeeze(DataMat_norm(1,:,:,3:end)).*repmat(mask,1,1,348);
plot(mean(reshape(DataMat_test1,75*120, 348), 1)), hold on

DataMat_test2 = squeeze(DataMat_norm(2,:,:,1:end-2)).*repmat(mask,1,1,348);
plot(mean(reshape(DataMat_test2,75*120, 348), 1))