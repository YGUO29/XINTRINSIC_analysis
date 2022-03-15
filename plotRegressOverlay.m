function plotRegressOverlay(output, para)

%% overlay normalized regression map for frequency, temp mod, spec mod in RGB channel
figurex([669,738,1263,420]); 
img_rgb = zeros(para.height, para.width, 3);
for i = 1:3
    b_ind = reshape(output.Rsquared(:,i),[para.height, para.width]);
    b_ind = b_ind-min(b_ind(:));
%     temp = temp./max(temp(:));
    img_rgb(:,:,i) = b_ind;
end
img_rgb = img_rgb./max(img_rgb(:));

subplot(1,2,1)
imagesc(img_rgb), axis image, axis off
title('R=freq, G=temp.mod., B=spec.mod.')
if isfield(para, 'ct')
    plotContour(para.ct);
end

if para.mirror
    set(gca,'XDir','reverse')
end


% overlay MAX of regression map for frequency, temp mod, spec mod in RGB channel
img_rgb = zeros(para.height, para.width, 3);
b_ind = reshape(output.Rsquared(:,1:3),[para.height, para.width, 3]);
[~, ind_max] = max(b_ind,[],3);
for i = 1:para.height
    for j = 1:para.width
%         img_rgb(i,j,ind_max(i,j)) = 1;
        img_rgb(i,j,ind_max(i,j)) = b_ind(i,j, ind_max(i,j));
    end
end
img_rgb = (img_rgb - min(img_rgb(:)))./range(img_rgb(:));
subplot(1,2,2)

imagesc(img_rgb), axis image, axis off
title('R=freq, G=temp.mod., B=spec.mod.')
if isfield(para, 'ct')
    plotContour(para.ct);
end

if para.mirror
    set(gca,'XDir','reverse')
end


