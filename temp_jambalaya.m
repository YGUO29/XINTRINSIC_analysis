% temporary for processing jambalaya session
filename = 'D:\=sounds=\Vocalization\Voc_jambalaya\Sound info.xlsx';
[Sd.wav, Sd.fs] = audioread('X:\=Sounds=\Vocalization\Sound_jambalaya_voc_(3.0+243.9+3.1)s_30dBatt.wav');
% 0.5s silence attached after each sound segment
opts = detectImportOptions(filename);
opts.SelectedVariableNames = 2; 
Duration = readmatrix(filename,opts);
opts.SelectedVariableNames = 1; 
N = readmatrix(filename,opts);
%% give each sound a number label
Names = cell(length(N),2);
% for i = 1:length(N)
%     temp = strsplit(N{i},' ');
%     Names{i,1} = lower(temp{4});
%     if contains(Names{i,1},'phee')
%         Names{i,2} = [Names{i,2}, 1];
%     end
%     if contains(Names{i,1},'tril')
%         Names{i,2} = [Names{i,2}, 2];
%     end
%     if contains(Names{i,1},'twit') 
%         Names{i,2} = [Names{i,2}, 3];
%     end
%     if contains(Names{i,1},'peep')
%         Names{i,2} = [Names{i,2}, 4];
%     end
%     if contains(Names{i,1},'tsik')
%         Names{i,2} = [Names{i,2}, 5];
%     end
% end
for i = 1:length(N)
    temp = strsplit(N{i},' ');
    temp2 = strsplit(temp{4},'_');
    Names{i,1} = lower(temp2{1});
    if strcmp(Names{i,1},'phee') | strcmp(Names{i,1},'pheestr') | strcmp(Names{i,1},'pheestrg')
        Names{i,2} = [Names{i,2}, 1];
    elseif strcmp(Names{i,1},'trill')
        Names{i,2} = [Names{i,2}, 2]; 
    elseif strcmp(Names{i,1},'trillphee') | strcmp(Names{i,1},'trilphee')
        Names{i,2} = [Names{i,2}, 3];
    elseif strcmp(Names{i,1},'twitter') | strcmp(Names{i,1},'twitter1')
        Names{i,2} = [Names{i,2}, 4];
    else 
        Names{i,2} = [Names{i,2}, 5];
    end
end
Names = Names(1:2:end, :);
% ind_single = cellfun(@(x) length(x)==1, Names(:,2), 'UniformOutput', false);
% ind_single = cell2mat(ind_single);

%% segregate response temporal trace according to stimuli
opt = struct;
opt.ampLimit    = 0.25.*[0 1];
[X, DataMat_norm] = getX(DataMat, para, opt);
% [X, DataMat_norm] = ViewData(DataMat, para, opt); 

DataMat_norm = squeeze(DataMat_norm);
tStart = para.preStim; % exclude pre-stimulus time
DataMat_sep = cell(1,length(N)/2);
Sd_sep = cell(1,length(N)/2);
for i = 1:length(N)/2
    ind = 2*(i-1)+1; % the ith stimulus, ind-th segment
    dur = 2*Duration(ind) + 1; % read this duration of data (with 1s silence) 
    nFrame = round(dur*para.fr);
    nSample = round(dur*Sd.fs);
    
    DataMat_sep{i} = squeeze(mean(DataMat_norm(:,:,round(tStart*para.fr)+1:round(tStart*para.fr)+nFrame), 3));
    % ===== normalize to [-1, 1] =====
    DataMat_sep{i} = DataMat_sep{i}./max(max(abs(DataMat_sep{i})));
    % ================================
    Sd_sep{i} = Sd.wav( round(tStart*Sd.fs)+1:round(tStart*Sd.fs)+nSample );
    tStart = tStart + dur;
end

%% Visualize response patterns seperately, generate data matrix "X"
X = zeros(length(DataMat_sep), size(DataMat_norm,1)*size(DataMat_norm,2));

DataMat_group = repmat(zeros(para.height, para.width), 1, 1, 5);
nDataMat = zeros(1,5);
figure, set(gcf,'Color', 'White', 'Position', [153         213        1806         858]);
[p,n] = numSubplots(size(Names,1));
ha = tight_subplot(p(1),p(2),[0 0],[.1 .01],[.01 .01]);
Max = max( cellfun(@(x) max(max(abs(x))), DataMat_sep) );

ind = 1:length(DataMat_sep);
% ind = outperm;
for i = 1:length(DataMat_sep)
    axes(ha(i));
    imagesc(DataMat_sep{ind(i)},[-Max Max]), axis image, axis off,colorbar off
    colormap(jet)
    title(Names{ind(i),1})
%     pause
    X(i,:) = reshape(DataMat_sep{ind(i)}, 1, size(DataMat_norm,1)*size(DataMat_norm,2));
    
    label = Names{ind(i),2};
    if ~isempty(label)    
        for k = 1:length(label)
            DataMat_group(:,:,label(k)) = DataMat_group(:,:,label(k)) + DataMat_sep{ind(i)};
            nDataMat(label(k)) = nDataMat(label(k)) + 1;
        end
    end
end
DataMat_group_avg = DataMat_group;


for i = 1:size(DataMat_group,3)
DataMat_group_avg(:,:,i) = squeeze(DataMat_group(:,:,i))./nDataMat(i);
end
getColorbar(Max,'percent')

%% plot grouped response patterns
% Title = {'phee', 'trill', 'twitter', 'peep', 'tsik'};
Title = {'phee', 'trill', 'trillphee', 'twitter', 'others'};

Max = max(abs(DataMat_group_avg(:)));
Min = min(abs(DataMat_group_avg(:)));
figure, set(gcf, 'Color', 'White');

for i = 1:5
    subplot(1,5,i), imagesc(DataMat_group_avg(:,:,i), [-Max, Max]), axis image
    colormap(jet)
    title(Title(i))
end
getColorbar(Max,'percent')
% for i = 1:5
%     mask_outline = zeros(size(DataMat_group_avg(:,:,i)));
%     mask = DataMat_group_avg(:,:,i)>0.7*max(DataMat_group_avg(:,:,i));
%     mask_outline = boundarymask(mask,4);
%     % figure, imshow(mask_outline);
%     figure, imshow(~mask_outline), axis image;
% end
%% subtractive analysis
mode = 'original'; % normalize or original
ind = nchoosek(1:5,2);
figure, set(gcf, 'Color', 'White', 'Position', [783         254        2068         564]);
[p,n] = numSubplots(size(ind,1));
ha = tight_subplot(p(1),p(2),[0.02 0.02],[.1 .01],[.01 .01]);
for i = 1:size(ind,1)
    switch mode
        case 'original'
            temp = DataMat_group_avg(:,:,ind(i,1)) - DataMat_group_avg(:,:,ind(i,2));
        case 'normalize'
            temp = DataMat_group_avg(:,:,ind(i,1))./max(max(abs(DataMat_group_avg(:,:,ind(i,1)))))...
                - DataMat_group_avg(:,:,ind(i,2))./max(max(abs(DataMat_group_avg(:,:,ind(i,2)))));
        otherwise
    end
    Max = max(abs(temp(:)));
    axes(ha(i));
    imagesc(temp, [-Max, Max]), axis image, axis off, colormap(jet), colorbar
    title([Title{ind(i,1)}, '-', Title{ind(i,2)}])
%     pause
end
getColorbar(Max,'number')
%% clustering of responses
% distance matrix method: euclidean
Dist = pdist(X);
Dist_mat = squareform(Dist);
Tree = linkage(Dist, 'average'); 
coph = cophenet(Tree, Dist)
Clus = cluster(Tree,'MaxClust',5);
[H, T, outperm] = dendrogram(Tree, 0);

% generate distance matrix, rearrange accoring to dendrogram
X_perm = X(outperm,:);
Dist_perm = pdist(X_perm);
Dist_perm_mat = squareform(Dist_perm);


figure, set(gcf, 'Color', 'White', 'Position', [618         364        1337         974]);
subplot(2,3,[1 4]), % dendrogram
    [H, T, outperm] = dendrogram(Tree, 0,'Orientation','Left','ColorThreshold','default');
    set(gca,'xtick',[])
    title('Dendrogram - response pattern')
ax1 = subplot(2,3,3); % response similarity   
    imagesc(flipud(Dist_perm_mat)), colormap(ax1, hsv), axis image, axis off
    title('Similarity matrix - response pattern')
ax2 = subplot(2,3,[2 5]); % response similarity with labels
    imagesc(flipud(Dist_perm_mat)), colormap(ax2, hsv), colorbar
    Names_perm = Names(outperm, :);
    set(gca,'ytick',1:66, 'yticklabels',flipud(Names_perm(:,1)));
ax3 = subplot(2,3,6); % label similarity
    labels = cell2mat(Names_perm(:,2));
    Dist_label = pdist(labels);
    Dist_label_mat = squareform(Dist_label);
    imagesc(flipud(logical(Dist_label_mat))), axis image, axis off, colormap(ax3, jet)
    title('Similarity matrix - vocalization type')
    
%% clustering of labels, and plot correponding reponse similarity matrix
% distance matrix method 2: euclidean
labels = cell2mat(Names(:,2));
Dist = pdist(labels);
Dist_mat = squareform(Dist);
Tree = linkage(Dist, 'single'); 
coph = cophenet(Tree, Dist)
Clus = cluster(Tree,'MaxClust',5);
[H, T, outperm] = dendrogram(Tree, 0);

% generate distance matrix, rearrange accoring to dendrogram
labels_perm = labels(outperm,:);
Dist_perm = pdist(labels_perm);
Dist_perm_mat = squareform(Dist_perm);


figure, set(gcf, 'Color', 'White', 'Position', [618         364        1337         974]);
subplot(2,3,[1 4]), % dendrogram
    [H, T, outperm] = dendrogram(Tree, 0,'Orientation','Left','ColorThreshold','default');
    set(gca,'xtick',[])
    title('Dendrogram - response pattern')
ax1 = subplot(2,3,3); % response similarity   
    imagesc(flipud(Dist_perm_mat)), colormap(ax1, hsv), axis image, axis off
    title('Similarity matrix - response pattern')
ax2 = subplot(2,3,[2 5]); % response similarity with labels
    imagesc(flipud(logical(Dist_perm_mat))), colormap(ax2, jet), colorbar
    Names_perm = Names(outperm, :);
    set(gca,'ytick',1:66, 'yticklabels',flipud(Names_perm(:,1)));
ax3 = subplot(2,3,6); % label similarity
    X_perm = X(outperm,:);
    Dist_label = pdist(X_perm);
    Dist_label_mat = squareform(Dist_label);
    imagesc(flipud(Dist_label_mat)), axis image, axis off, colormap(ax3, hsv)
    title('Similarity matrix - vocalization type')
%% load sound, display spectrograms in rearranged sequence
Spec = cell(1,length(Sd_sep));
ind = 1:length(Sd_sep);
for i = 1:length(Sd_sep)
    Sd_temp.fs = Sd.fs;
    Sd_temp.wav = Sd_sep{ind(i)};
%     axes(ha(i)); getSpectrogram(Sd,1,0.01);
    [~,Spec{i},~,~,~] = getCochleogram(Sd_temp, 0.01, 'ERB', 0);
%     rectangle('Position',[1, 1, size(Spec{i},2), size(Spec{i}, 1)],'EdgeColor','y');
    axis off
    colorbar off
end

%% plot cochleograms
L_pad = max(cell2mat( cellfun(@(x) size(x,2), Spec,'UniformOutput',false) ));
figure, set(gcf,'Color', 'White', 'Position', [153         213        1806         858]);
[p,n] = numSubplots(length(Sd_sep));
ha = tight_subplot(p(1),p(2),[.02 .005],[.1 .1],[.01 .01]);
ind = outperm;
for i = 1:length(Sd_sep)
    axes(ha(i)); 
    temp = [Spec{ind(i)}, zeros(size(Spec{ind(i)},1), L_pad - size(Spec{ind(i)},2))];
    imagesc(flipud(temp));
    rectangle('Position',[1, 1, size(Spec{ind(i)},2), size(Spec{ind(i)}, 1)],'EdgeColor','y');

    title(Names_perm{i})
    axis off
    colorbar off
    colormap(jet)
end
