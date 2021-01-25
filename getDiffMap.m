% compare two conditions - get differential maps and statistics
%% compare two sets of trials
% X1 = X(1:size(X,1)/2, :); X2 = X(size(X,1)/2+1:end, :);
X2 = X(2:2:24,:);
X1 = X(1:2:24,:); 

figurex;
% for i = 1:4
%     subplot(2,2,i)
%     eval(['temp = mean(X',num2str(i),', 1);']);
    img = mean(X1 - X2, 1); % X1: Original, X2: modified
    Max = max(abs(img));
    imagesc(reshape(img, para.height, para.width), [-Max, Max].*0.8), 
    axis image, colorbar
    colormap(jet)
%     title('Red: Trill-phee > Trill')
    % colormap(cbrewer('div', 'RdBu', 256))
% end
%% get statistics
[~,p,~,stats] = ttest2(X1, X2);

%% combine t-value and p-value
cmap = flipud(cbrewer('div', 'RdBu', 256));

figurex;
img = cell(1, 3);

% figure1: t-value only
img{1} = reshape(stats.tstat, para.height, para.width);

% figure2: p-value only
img_amp = -log10(p);
img{2} = reshape(img_amp, para.height, para.width);

% figure3: t-value and p-value together
cmap_ind = rescale(stats.tstat, 1, 256);
cmap_ind(ismissing(cmap_ind)) = 256/2;
img_rgb = cmap(floor(cmap_ind),:); % RGB values
img_amp(ismissing(img_amp)) = 0;
img_amp = rescale(img_amp);

img{3} = img_rgb.*repmat(img_amp', 1, 3);
img{3} = reshape(img{3}, para.height, para.width, 3);

Titles = {'t-values', '-log10(p-value)', 't-values with p-values as amplitude mask'};
for i = 1:3
    subplot(1,3,i)
    Max = max(abs(img{i}(:)));
    imagesc(img{i}, [-Max, Max].*0.8), 
    axis image, colorbar
    colormap(cmap)
    title(Titles{i})
end
