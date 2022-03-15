function plotMap_with2Dcolorbar(sat, hue, para)
% hue: 2d map of tuning
% sat: 2d map of saturation and value
% para.hue_max: max value of hue
% para.continuous_cb: 1 for continuous colorbar

map = repmat(hue, 1, 1, 3);
map(:,:,1) = (map(:,:,1)./para.hue_max).*0.8; % H
map(:,:,2) = sat; % S
map(:,:,2) = map(:,:,2)./max(max(map(:,:,2))); % S
map(:,:,3) = map(:,:,2); % V

map = hsv2rgb(map);
if para.mirror
    map = fliplr(map);
end

imagesc(map), axis image, axis off

cmap_hsv = ones(para.hue_max, 3);
if ~para.continuous_cb % discrete colorbar
    cmap_hsv(:,1) = ((1:para.hue_max)./para.hue_max).*0.8';
else
    cmap_hsv(:,1) = ((1:256)./256).*0.8';
end
cmap = hsv2rgb(cmap_hsv);

colormap(gca, cmap);
hcb = colorbar;
cb_ticks = range(hcb.Limits)/para.hue_max : ...
    range(hcb.Limits)/para.hue_max : hcb.Limits(2);
set(hcb, 'Ticks', cb_ticks, 'TickLabels',para.cb_labels);
if isfield(para, 'title')
    title(para.title)
end
end