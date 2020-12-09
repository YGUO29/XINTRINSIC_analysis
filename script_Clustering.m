% script for clustering
% X: 1 row = 1 observation

%% Select part of X to do clustering
X_orig = X; 
X = X_orig(sort(category_ind, 'ascend') + 15, :);

%% 
% distance matrix method 1: 1-corr
ind = nchoosek(1:size(X,1),2);
r = zeros(1,length(ind));
for i = 1:size(ind,1)
    rr = corrcoef(X(ind(i,1),:), X(ind(i,2),:));
    r(i) = rr(1,2);
end
Dist = 1-r;
%%
% distance matrix method 2: euclidean
Dist = pdist(X); % default: euclidean

Dist_mat = squareform(Dist);
Tree = linkage(Dist, 'average'); 
coph = cophenet(Tree, Dist)

% plot dendrogram
figure, [H, T, outperm] = dendrogram(Tree, 0);

% generate distance matrix, rearrange accoring to dendrogram
X_perm = X(outperm,:);
Dist_perm = pdist(X_perm);
Dist_perm_mat = squareform(Dist_perm);
figure, imagesc(Dist_perm_mat), colormap(jet)

%% display data in lower dimension
nDim = 2;
Y = mdscale(Dist_mat,nDim);
Clus = cluster(Tree,'MaxClust',5); 

figure, 
switch nDim
    case 2
    scatter(Y(:,1), Y(:,2), 36, Clus,'filled');
    case 3
    scatter3(Y(:,1), Y(:,2), Y(:,3), 36, Clus, 'filled');
end
%% DISPLAY
figure
cmax = 0.08;
cmin = -0.08;

for i = 1:length(outperm)
    subplot(4,13,i)
    imagesc(reshape(-X_perm(i,:), para.height, para.width), [cmin, cmax])
%     imagesc(reshape(-X(i,:), para.height, para.width), [cmin, cmax])
    axis off
    axis image
    colormap(jet)
end
%% DISPLAY responses
% dynamic data - in the order or rearranged sequence (similar responses are close to each other)
tic
opt = struct;
opt.ampLimit    = 0.08.*[-1, 1];
opt.trials = outperm; 
figure, set(gcf, 'color','w')
[X, DataMat_norm] = ViewData(DataMat, para, opt); % X may contain NaNs if there are masked pixels
%     [X, DataMat_norm] = getX(DataMat, para, opt);
toc

% ======== construct a X without NaN ========
[~, ind_delete]     = find( isnan(X) ); % linear index
ind_save            = setdiff(1:para.width*para.height, ind_delete);
X(:, ind_delete)    = [];

%% load sound, display spectrograms in rearranged sequence
soundpath = 'D:\SynologyDrive\=sounds=\Natural sound\Natural_JM_ModelMatched_Sounds\MusicSpeech\';
addpath(genpath(soundpath))
list_full = dir(fullfile(soundpath,'*.wav'));
names = natsortfiles({list_full.name})';
names_new = names(outperm);

[p,n] = numSubplots(length(names_new));
ha = tight_subplot(p(1),p(2),[0 0],[.1 .01],[.01 .01]);
figure,
Spec = zeros(202, 200, length(names_new));
for i = 1:length(names_new)
    [Sd.wav, Sd.fs] = audioread(names_new{i});
%     axes(ha(i)); getSpectrogram(Sd,1,0.01);
    axes(ha(i)); [~,Spec(:,:,i),~,~,~] = getCochleogram(Sd, 0.01, 'ERB', 0);
    axis off
    colorbar off
end

%% distance matrix method 2: euclidean
Dist = pdist(X_spec);
Dist_mat = squareform(Dist);
Tree = linkage(Dist, 'centroid'); 
coph = cophenet(Tree, Dist)
figure, [H, T, outperm] = dendrogram(Tree, 0);
Clus = cluster(Tree,'MaxClust',8);
% generate distance matrix, rearrange accoring to dendrogram
X_spec_perm = X_spec(outperm,:);
Dist_perm = pdist(X_spec_perm);
Dist_perm_mat = squareform(Dist_perm);
figure, imagesc(Dist_perm_mat), colormap(jet)
%%
load('D:\=code=\Sound_analysis\F_yg_marm.mat') 
X_feature = F.F_mat';
X_feature_perm = X_feature(outperm,:);
Dist_perm = pdist(X_feature_perm);
Dist_perm_mat = squareform(Dist_perm);
figure, imagesc(Dist_perm_mat), colormap(jet)
%%
load('D:\SynologyDrive\=data=\category_regressors_withLZvoc.mat')
C = C_voc;
tags = C.category_assignments; 
tags_inorder = tags(outperm);
nTags = max(tags);
Color = C.colors;
Color_inorder = Color()
figure, 
b = bar(1:nTags, ones(1,nTags),'FaceColor','flat')
b.CData = Color_inorder;

