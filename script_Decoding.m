% decoding for natural sounds
% load('D:\SynologyDrive\=data=\category_regressors.mat')
% folder_sound = 'D:\SynologyDrive\=code=\McdermottLab\sound_natural\';
load('D:\SynologyDrive\=data=\category_regressors_LZVoc_major.mat')
folder_sound = 'D:\SynologyDrive\=sounds=\Vocalization\temp_for capsule\AllMajor_Orig_norm';
list = dir(fullfile(folder_sound,'*.wav'));
[~, ind, ~]     = natsortfiles({list.name});
list = list(ind);
%% Load data (DataMat)
% animal = '96B'; 
% session = 'NatVoc';
% datapath = 'D:\SynologyDrive\=data=\XINTRINSIC';
% load([datapath, '\', animal, '\DataMat_', animal, '_', session, '_', num2str(para.nRep), 'reps.mat'])
% [DataMat, para] = getDataMat;

[file, datapath] = uigetfile('D:\SynologyDrive\=data=\XINTRINSIC\*.mat');
load(fullfile(datapath, file))

opt             = struct;
% opt.tWindow     = [para.preStim, para.preStim + para.durStim + 4];
opt.tWindow     = [para.preStim.*ones(para.nStim,1), para.preStim.*ones(para.nStim,1) + dur_mat'];
[X1, ~, mov_rel_sep] = getX(DataMat, para, opt);
% X_sep: [nRep, nStim, height, width, frames]
% X_sep = X_sep(:,:,:,20:110, :); % cut part of the image (width)
para.width1 = size(mov_rel_sep,4);

% include temporal info
% X_sep = X_sep(:,:,:,:,para.fr*para.preStim+1:3:para.fr*(para.preStim+para.durStim));
% X_sep = reshape(X_sep, size(X_sep,1),  size(X_sep,2), size(X_sep,3)*size(X_sep,4)* size(X_sep,5));

% don't include temporal info
X_sep = zeros(para.nRep, para.nStim, para.height*para.width);
for i = 1:para.nStim
temp = squeeze(mean( mov_rel_sep(:,i,:,:,floor(para.fr*opt.tWindow(i,1))+1 : floor(para.fr*opt.tWindow(i,2))) , 5 ));
X_sep(:,i,:) = reshape(temp, size(temp,1),  size(temp,2)*size(temp,3));
end
% X_sep = reshape(X_sep, size(X_sep,1)*size(X_sep,2),  size(X_sep,3));
% X = X_sep;

%% apply mask
mask = reshape(mask1(:,20:110),75*91,1);
mask = repmat(mask,[1,10, 165]);
mask = permute(mask,[2 3 1]);
X = X_sep.*~mask;

%% start decoding procedure
para_decoding.variable = 'name'; % name, cat, or sub
para_decoding.shuffle = 0;

tic
for i = 1:para.nRep % i: the repetition for TEST
    % all other repetitions for TRAIN
    ind_tr = setdiff(1:para.nRep,i);
    Tr = table;
    Tr.resp     = reshape(X_sep(ind_tr,:,:), (para.nRep-1)*size(X_sep,2), size(X_sep,3));
    Tr.number   = repelem(1:para.nStim, para.nRep-1)';
    
    % names
    nametemp    = {list.name};
    if para_decoding.shuffle
        nametemp = nametemp(randperm(length(nametemp))); % shuffle labels
    end      
    Tr.name      = repelem(nametemp, para.nRep-1)';
    
    % category labels       
    if size(C.category_labels,1) < size(C.category_labels,2)
        C.category_labels = C.category_labels';
    end
    Tr.category  = repelem(C.category_labels(C.category_assignments), para.nRep-1);
    if para_decoding.shuffle
        Tr.category  = Tr.category(randperm(length(Tr.category))); % shuffle labels
    end
    
    % subject labels
    if size(C.subject_labels,1) < size(C.subject_labels,2)
        C.subject_labels = C.subject_labels';
    end
    Tr.subject  = repelem(C.subject_labels(C.subject_assignments), para.nRep-1);
    if para_decoding.shuffle
        Tr.subject  = Tr.subject(randperm(length(Tr.subject))); % shuffle labels
    end
    
    switch para_decoding.variable
        case 'name'
        % train on identification task
        md = fitcecoc(Tr.resp,Tr.name); % adding NumNeighbors did not improve results
%         md = fitcknn(Tr.resp, Tr.name);
        
         % test dataset
        Te = table;
        Te.resp = squeeze(X_sep(i,:,:));
        Te.predname = predict(md, Te.resp);
        Te.truename = natsortfiles({list.name})';
        PercCorrect(i) = sum(cellfun(@strcmp, Te.predname, Te.truename))./size(Te,1);
        
        case 'cat'
        % train on classification task
        md = fitcecoc(Tr.resp,Tr.category); % adding NumNeighbors did not improve results
%         md = fitcknn(Tr.resp, Tr.category);
        Te = table;
        Te.resp = squeeze(X_sep(i,:,:));
        Te.predcategory = predict(md, Te.resp);
        Te.truecategory = C.category_labels(C.category_assignments);
        PercCorrect(i) = sum(cellfun(@strcmp, Te.predcategory, Te.truecategory))./size(Te,1);

        case 'sub'
        md = fitcecoc(Tr.resp, Tr.subject);
%         md = fitcknn(Tr.resp, Tr.subject);
        Te = table;
        Te.resp = squeeze(X_sep(i,:,:));
        Te.predsubject = predict(md, Te.resp);
        Te.truesubject = C.subject_labels(C.subject_assignments);
        PercCorrect(i) = sum(cellfun(@strcmp, Te.predsubject, Te.truesubject))./size(Te,1);
   
        otherwise
    end

%     cm = confusionchart(Te.truename, Te.predname,'RowSummary','row-normalized');
%     ind_linear = find(cm.NormalizedValues);
%     [ii,jj] = ind2sub(size(cm.NormalizedValues), ind_linear);
%     PercCorrect(2,i) = length(find(ii-jj == 0))./size(cm.NormalizedValues,1); 
    i
end
toc
PercCorrect_all(2,:) = PercCorrect;
% 1,2: shuffled and name and shuffled
% 3,4: shuffled and cat
%%
figurex([335         192        1364         654]); hold on
boxplot(PercCorrect_all', 'Colors','bbrrkk')
line([0 8], [1/para.nStim, 1/para.nStim], 'LineStyle', '--')
line([0 8], [1/length(C.category_labels), 1/length(C.category_labels)],...
    'LineStyle', '--', 'color', 'r')
line([0 8], [1/length(C.subject_labels), 1/length(C.subject_labels)],...
    'LineStyle', '--', 'color', 'k')

set(findobj(gca,'type','line'),'linew',2)
% xticklabels({'Nearest neighbor', 'SVM', 'Nearest neighbor', 'SVM'})
% xticklabels({'all(NN)', 'core(NN)', 'non-core(NN)', 'all(SVM)','core(SVM)','non-core(SVM)', 'Shuffled'})
xticklabels({'shuffle-sound', 'sound', 'shuffle-category', 'category','shuffle-subject','subject'})

ylabel('correct rate')
ylim([0 1])
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