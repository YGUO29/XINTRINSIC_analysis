function plotContour(ct, color)
% ct is the contour object
% set color to be 'jet', 'RdBu' or 'gray'

% figurex; 
hold on
% [~, ~, ic] = unique(ct.level);

switch color
    case 'jet'
        cmap = jet(length(ct.zindex));
    case 'RdBu'
        cmap = cbrewer('div', 'RdBu', length(ct.zindex));
        cmap = flipud(cmap);
    case 'gray'
        cmap = gray(length(ct.zindex));
    otherwise
end


for i = 1:length(ct.level)
    [~, cmap_ind] = min(abs(ct.level(i) - ct.zindex));
    plot(ct.c{i}(1,:), ct.c{i}(2,:), 'color', cmap(cmap_ind,:),'linewidth',1.5)
%     plot(ct.c{i}(1,:), ct.c{i}(2,:), 'color', cmap(cmap_ind,:))
end

end