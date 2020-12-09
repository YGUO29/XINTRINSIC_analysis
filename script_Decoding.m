% decoding for natural sounds
% load('D:\SynologyDrive\=data=\category_regressors.mat')
% folder_sound = 'D:\SynologyDrive\=code=\McdermottLab\sound_natural\';
load('D:\SynologyDrive\=data=\category_regressors_withLZvoc.mat')
folder_sound = 'D:\SynologyDrive\=sounds=\Natural sound\Natural_JM_XINTRINSIC_withLZVoc_200909\Norm';
list = dir(fullfile(folder_sound,'*.wav'));

%% Load data (DataMat)
animal = '96B'; 
session = 'NatVoc';
datapath = 'D:\SynologyDrive\=data=\XINTRINSIC';
load([datapath, '\', animal, '\DataMat_', animal, '_', session, '_', num2str(para.nRep), 'reps.mat'])

% [DataMat, para] = getDataMat;

opt             = struct;
opt.tWindow     = [para.preStim, para.preStim + para.durStim + 4];
[X1, ~, X_sep] = getX(DataMat, para, opt);
% X_sep: [nRep, nStim, height, width, frames]
% X_sep = X_sep(:,:,:,20:110, :); % cut part of the image (width)
para.width1 = size(X_sep,4);

% include temporal info
% X_sep = X_sep(:,:,:,:,para.fr*para.preStim+1:3:para.fr*(para.preStim+para.durStim));
% X_sep = reshape(X_sep, size(X_sep,1),  size(X_sep,2), size(X_sep,3)*size(X_sep,4)* size(X_sep,5));

% don't include temporal info
X_sep = squeeze(mean( X_sep(:,:,:,:,floor(para.fr*opt.tWindow(1))+1 : floor(para.fr*opt.tWindow(2))) , 5 ));
X_sep = reshape(X_sep, size(X_sep,1),  size(X_sep,2),  size(X_sep,3)* size(X_sep,4));
% X_sep = reshape(X_sep, size(X_sep,1)*size(X_sep,2),  size(X_sep,3));

%% apply mask
mask = reshape(mask1(:,20:110),75*91,1);
mask = repmat(mask,[1,10, 165]);
mask = permute(mask,[2 3 1]);
X = X_sep.*~mask;

%% create a table for training data 
PercCorrect2 = zeros(1,10);
%% 
tic
for i = 1:para.nRep % i: the repetition for TEST
    % all other repetitions for TRAIN
    ind_tr = setdiff(1:para.nRep,i);
    Tr = table;
    Tr.resp      = reshape(X(ind_tr,:,:), (para.nRep-1)*size(X,2), size(X,3));
    Tr.number    = repelem(1:para.nStim, para.nRep-1)';
    nametemp     = natsortfiles({list.name});
%     nametemp     = nametemp(randperm(length(nametemp))); % shuffle labels
    Tr.name      = repelem(nametemp,para.nRep-1)';
    Tr.category  = repelem(C.category_labels(C.category_assignments), para.nRep-1); 
%     Tr.category  = Tr.category(randperm(length(Tr.category)));
% Tr.colors    = C.colors(C.category_assignments,:);

    % train on identification task
%     md = fitcecoc(Tr.resp,Tr.name, 'Learners', 'knn'); % adding NumNeighbors did not improve results
    md = fitcknn(Tr.resp,Tr.name);
    % train on classification task
%     md = fitcecoc(Tr.resp,Tr.category); % adding NumNeighbors did not improve results
%     md = fitcknn(Tr.resp,Tr.category);
    
    Te = table;
    Te.resp = squeeze(X(i,:,:));
    Te.predname = predict(md,Te.resp);
    Te.truename = natsortfiles({list.name})';
    PercCorrect1(4,i) = sum(cellfun(@strcmp, Te.predname, Te.truename))./size(Te,1);
%     Te.predcategory = predict(md,Te.resp);
%     Te.truecategory = C.category_labels(C.category_assignments);
%     PercCorrect2(1,i) = sum(cellfun(@strcmp, Te.predcategory, Te.truecategory))./size(Te,1);

%     cm = confusionchart(Te.truename, Te.predname,'RowSummary','row-normalized');
%     ind_linear = find(cm.NormalizedValues);
%     [ii,jj] = ind2sub(size(cm.NormalizedValues), ind_linear);
%     PercCorrect(2,i) = length(find(ii-jj == 0))./size(cm.NormalizedValues,1); 
    i
end
toc
%%
figurex; hold on
boxplot(PercCorrect1')
% line([0 8], [1/165 1/165], 'LineStyle', '--')
line([0 8], [1/11 1/11], 'LineStyle', '--')

set(findobj(gca,'type','line'),'linew',2)
% xticklabels({'Nearest neighbor', 'SVM', 'Nearest neighbor', 'SVM'})
xticklabels({'all(NN)', 'core(NN)', 'non-core(NN)', 'all(SVM)','core(SVM)','non-core(SVM)', 'Shuffled'})

ylabel('correct rate')
ylim([0 0.7])
%% create a model
% md = fitcecoc(Tr.resp,Tr.name);

function Tr = ReduceDimension(X_sep)
X = reshape(X_sep,size(X_sep,1)*size(X_sep,2), size(X_sep,3));
[U, S, V] = svd(X);
eigval = diag(S);
figurex,
plot(cumsum(eigval(2:end))./sum(eigval(2:end)))
N = 1000;
X = U*S*V(:,1:N);
X_sep = reshape(X, size(X_sep,1), size(X_sep,2), size(X,2));

figurex,
for i = 1:6
    subplot(2,3,i)
    imagesc(reshape(V(:,i),para.height, para.width1));
    axis image, colorbar
end
end