function plotContour(ct)

% figurex; 
hold on
[~, ~, ic] = unique(ct.level);
% COLOR = hsv(length(ct.zindex));
COLOR = gray(length(ct.zindex));
COLOR = flipud(COLOR);
% COLOR = cbrewer('seq', 'YlOrRd', length(ct.zindex));
% 
for i = 1:length(ct.level)
    plot(ct.c{i}(1,:), ct.c{i}(2,:), 'color', COLOR(ic(i),:))
end

end