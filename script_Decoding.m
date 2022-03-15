% decoding for natural sounds
load('D:\SynologyDrive\=data=\F_halfcosine_marm_NatJM.mat')
% folder_sound = 'D:\SynologyDrive\=code=\McdermottLab\sound_natural\';
% load('D:\SynologyDrive\=data=\category_regressors_LZVoc_major.mat')
if ~isfield(F,'sound_names')
%     folder_sound = 'D:\SynologyDrive\=sounds=\Natural sound\Natural_JM_original';
    folder_sound = 'D:\SynologyDrive\=sounds=\Vocalization\temp_for capsule\AllMajor_Orig_norm2';
    list = dir(fullfile(folder_sound,'*.wav'));
    [~, ind, ~]     = natsortfiles({list.name});
    list = list(ind);
else
end

%% get reduced dimension of X_sep
R_sep = zeros(para.nRep, para.nStim, 6);
for i = 1:para.nRep
R_sep(i,:,:) = squeeze(X_sep(i,:,:)) * pinv(Decomp_102D.Ws{1});
end
%% Load data (DataMat)
% DataMat = DataMat(:,:,11:140,31:210,:);
% [DataMat, para] = getDataMat;

% [file, datapath] = uigetfile('D:\SynologyDrive\=data=\XINTRINSIC\*.mat');
% load(fullfile(datapath, file))

opt = struct;
opt.tWindow = [para.preStim, para.preStim + para.durStim];
% opt.tWindow = repmat(opt.tWindow, para.nStim, 1); % tWindow could be different for different trials
opt.tWindow = [para.preStim.*ones(para.nStim,1), para.preStim.*ones(para.nStim,1) + dur_mat'];
[X1, ~, mov_rel_sep] = getX(DataMat, para, opt);
% X_sep: [nRep, nStim, height, width, frames]
% X_sep = X_sep(:,:,:,20:110, :); % cut part of the image (width)
% para.width1 = size(mov_rel_sep,4);

% don't include temporal info
% conver to a matrix of #reps x #trials x #pixels
X_sep = zeros(para.nRep, para.nStim, para.height*para.width);
for i = 1:para.nStim
temp = squeeze(mean( mov_rel_sep(:,i,:,:,floor(para.fr*opt.tWindow(i,1))+1 : floor(para.fr*opt.tWindow(i,2))) , 5 ));
X_sep(:,i,:) = reshape(temp, size(temp,1),  size(temp,2)*size(temp,3));
end


% include temporal info
% X_sep = X_sep(:,:,:,:,para.fr*para.preStim+1:3:para.fr*(para.preStim+para.durStim));
% X_sep = reshape(X_sep, size(X_sep,1),  size(X_sep,2), size(X_sep,3)*size(X_sep,4)* size(X_sep,5));


%% apply mask
mask = reshape(mask1(:,20:110),75*91,1);
mask = repmat(mask,[1,10, 165]);
mask = permute(mask,[2 3 1]);
X = X_sep.*~mask;

%% start decoding procedure
X_input = R_sep;

para_decoding.variable = 'name'; % name, cat, or sub
para_decoding.shuffle = 0;

tic
for i = 1:para.nRep % i: the repetition for TEST
    % set up Train dataset: all other repetitions for TRAIN
    ind_tr = setdiff(1:para.nRep,i);
    Tr = table;
    Tr.resp = reshape(X_input(ind_tr,:,:), (para.nRep-1)*size(X_input,2), size(X_input,3));
    Tr.number = repelem(1:para.nStim, para.nRep-1)';
    
    % set up test dataset
    Te = table;
    Te.resp = squeeze(X_input(i,:,:));
    
    switch para_decoding.variable
        case 'name'
        % train on identification task
        % names
        if ~isfield(F,'sound_names')
            F.sound_names = {list.name};
        end
        name_temp = F.sound_names;

        if para_decoding.shuffle
            name_temp = name_temp(randperm(length(name_temp))); % shuffle labels
        end

        if size(name_temp,2)>size(name_temp,1)
            Tr.name = repelem(name_temp, para.nRep-1)';
        else
            Tr.name = repelem(name_temp, para.nRep-1);
        end
        Te.truename = F.sound_names';
        for j = 1:4 % 5 classifiers
            switch j
                case 5 
                % ==== LDA ====
%                 md = fitcdiscr(Tr.resp, Tr.name, 'discrimType', 'pseudoLinear'); % error showed as "Predictor x21 has zero within-class variance", so added pseudoLinear option
%                 Te.predname = predict(md, Te.resp); 
                case 1
                % ==== SVM ====
                md = fitcecoc(Tr.resp, Tr.name); 
                Te.predname = predict(md, Te.resp); 
                case 2
                % ==== 1NN ====
                md = fitcknn(Tr.resp, Tr.name, 'NumNeighbors',1); 
                Te.predname = predict(md, Te.resp); 
                case 3
                % ==== DT ====
                md = fitctree(Tr.resp, Tr.name); 
                Te.predname = predict(md, Te.resp); 
                case 4
                % ==== RF ====
                md = TreeBagger(50, Tr.resp, Tr.name, 'Method', 'classification', 'OOBPrediction','On');
                Te.predname = predict(md, Te.resp); 
                otherwise
            end

            PercCorrect(i,j) = sum(cellfun(@strcmp, Te.predname, Te.truename))./size(Te,1);
            [i, j]
        end
%         md = fitcecoc(Tr.resp,Tr.name); % adding NumNeighbors did not improve results

        
        case 'cat'
        % train on classification task
            % category labels       
        if isfield(F,'C'); C = F.C; end
        if size(C.category_labels,1) < size(C.category_labels,2)
            C.category_labels = C.category_labels';
        end
        category_temp = C.category_labels(C.category_assignments)';
        if para_decoding.shuffle
            category_temp  = category_temp(randperm(length(category_temp))); % shuffle labels
        end
        Tr.category  = repelem(category_temp, para.nRep-1)';

        % subject labels
        if size(C.subject_labels,1) < size(C.subject_labels,2)
            C.subject_labels = C.subject_labels';
        end
        Tr.subject  = repelem(C.subject_labels(C.subject_assignments), para.nRep-1);
        if para_decoding.shuffle
            Tr.subject  = Tr.subject(randperm(length(Tr.subject))); % shuffle labels
        end
        
        Te.truecategory = C.category_labels(C.category_assignments);
        for j = 1:4 % 5 classifiers
            switch j
                case 5 
                % ==== LDA ====
%                 md = fitcdiscr(Tr.resp, Tr.name, 'discrimType', 'pseudoLinear'); % error showed as "Predictor x21 has zero within-class variance", so added pseudoLinear option
%                 Te.predname = predict(md, Te.resp); 
                case 1
                % ==== SVM ====
                md = fitcecoc(Tr.resp, Tr.category); 
                Te.predcategory = predict(md, Te.resp); 
                case 2
                % ==== 1NN ====
                md = fitcknn(Tr.resp, Tr.category, 'NumNeighbors',1); 
                Te.predcategory = predict(md, Te.resp); 
                case 3
                % ==== DT ====
                md = fitctree(Tr.resp, Tr.category); 
                Te.predcategory = predict(md, Te.resp); 
                case 4
                % ==== RF ====
                md = TreeBagger(50, Tr.resp, Tr.category, 'Method', 'classification', 'OOBPrediction','On');
                Te.predcategory = predict(md, Te.resp); 
                otherwise
            end

            PercCorrect(i,j) = sum(cellfun(@strcmp, Te.predcategory, Te.truecategory))./size(Te,1);
            [i, j]
        end
        
        case 'sub'
        Te.truesubject = C.subject_labels(C.subject_assignments);

        for j = 1:4 % 5 classifiers
            switch j
                case 5 
                % ==== LDA ====
%                 md = fitcdiscr(Tr.resp, Tr.name, 'discrimType', 'pseudoLinear'); % error showed as "Predictor x21 has zero within-class variance", so added pseudoLinear option
%                 Te.predname = predict(md, Te.resp); 
                case 1
                % ==== SVM ====
                md = fitcecoc(Tr.resp, Tr.subject); 
                Te.predsubject = predict(md, Te.resp); 
                case 2
                % ==== 1NN ====
                md = fitcknn(Tr.resp, Tr.subject, 'NumNeighbors',1); 
                Te.predsubject = predict(md, Te.resp); 
                case 3
                % ==== DT ====
                md = fitctree(Tr.resp, Tr.subject); 
                Te.predsubject = predict(md, Te.resp); 
                case 4
                % ==== RF ====
                md = TreeBagger(50, Tr.resp, Tr.subject, 'Method', 'classification', 'OOBPrediction','On');
                Te.predsubject = predict(md, Te.resp); 
                otherwise
            end

            PercCorrect(i,j) = sum(cellfun(@strcmp, Te.predsubject, Te.truesubject))./size(Te,1);
            [i, j]
        end
        otherwise
            
    end

%     cm = confusionchart(Te.truename, Te.predname,'RowSummary','row-normalized');
%     ind_linear = find(cm.NormalizedValues);
%     [ii,jj] = ind2sub(size(cm.NormalizedValues), ind_linear);
%     PercCorrect(2,i) = length(find(ii-jj == 0))./size(cm.NormalizedValues,1); 
    
end
toc
% PercCorrect_all(1,:) = PercCorrect;
% 1,2: shuffled and name and shuffled
% 3,4: shuffled and cat

%% save data
decode_NN = table; decode_NN.sound = PercCorrect(:,2);
decode_SVM = table; decode_SVM.sound = PercCorrect(:,1);
decode_RF = table; decode_RF.sound = PercCorrect(:,4);
decode_DT = table; decode_DT.sound = PercCorrect(:,3);


% save('E:\Dropbox\_Yueqi Guo\_research\Imaging\=data=\102D\Decoding_Voc_2103014.mat',...
%     'decode_NN', 'decode_SVM', 'decode_RF', 'decode_DT');

figurex;
violinplot([decode_NN.sound, decode_SVM.sound, ...
    decode_DT.sound, decode_RF.sound]);
xticklabels({'1NN', 'SVM', 'DT', 'RF'})
ylabel('classification accuracy')
title('decoding the specific sound')
hold on,
line([0 4.5], [1/para.nStim, 1/para.nStim], 'LineStyle', '--')

decode_NN.category = PercCorrect(:,2);
decode_SVM.category = PercCorrect(:,1);
decode_RF.category = PercCorrect(:,4);
decode_DT.category = PercCorrect(:,3);
figurex;
violinplot([decode_NN.category, decode_SVM.category, ...
    decode_DT.category, decode_RF.category]);
xticklabels({'1NN', 'SVM', 'DT', 'RF'})
ylabel('classification accuracy')
title('decoding sound category')
hold on,
line([0 4.5], [1/4, 1/4], 'LineStyle', '--')

decode_NN.subject = PercCorrect(:,2);
decode_SVM.subject = PercCorrect(:,1);
decode_RF.subject = PercCorrect(:,4);
decode_DT.subject = PercCorrect(:,3);
figurex;
violinplot([decode_NN.subject, decode_SVM.subject, ...
    decode_DT.subject, decode_RF.subject]);
xticklabels({'1NN', 'SVM', 'DT', 'RF'})
ylabel('classification accuracy')
title('decoding the caller')
hold on,
line([0 4.5], [1/6, 1/6], 'LineStyle', '--')
%%
figurex([335         192        1364         654]); hold on
% boxplot(PercCorrect_NN')
boxplot(PercCorrect_SVM', 'Colors','bbrrkk')
line([0 8], [1/para.nStim, 1/para.nStim], 'LineStyle', '--')
line([0 8], [1/length(C.category_labels), 1/length(C.category_labels)],...
    'LineStyle', '--', 'color', 'r')
% line([0 8], [1/length(C.subject_labels), 1/length(C.subject_labels)],...
%     'LineStyle', '--', 'color', 'k')

set(findobj(gca,'type','line'),'linew',2)
xticklabels({'Sound (shuffled)', 'Sound', 'Category (shuffled)', 'Category'})
% xticklabels({'all(NN)', 'core(NN)', 'non-core(NN)', 'all(SVM)','core(SVM)','non-core(SVM)', 'Shuffled'})
% xticklabels({'shuffle-sound', 'sound', 'shuffle-category', 'category','shuffle-subject','subject'})
% xticklabels({'shuffle-sound', 'sound', 'shuffle-category', 'category'});
ylabel('correct rate')
ylim([0 1])
title('Decoding with SVM')

%%
figurex([335         192        1364         654]); hold on
% boxplot(PercCorrect_NN')
boxplot(PercCorrect_SVM', 'Colors','bbrrkk')
line([0 8], [1/para.nStim, 1/para.nStim], 'LineStyle', '--')
line([0 8], [1/length(C.category_labels), 1/length(C.category_labels)],...
    'LineStyle', '--', 'color', 'r')
% line([0 8], [1/length(C.subject_labels), 1/length(C.subject_labels)],...
%     'LineStyle', '--', 'color', 'k')

set(findobj(gca,'type','line'),'linew',2)
xticklabels({'Sound (shuffled)', 'Sound', 'Category (shuffled)', 'Category'})
% xticklabels({'all(NN)', 'core(NN)', 'non-core(NN)', 'all(SVM)','core(SVM)','non-core(SVM)', 'Shuffled'})
% xticklabels({'shuffle-sound', 'sound', 'shuffle-category', 'category','shuffle-subject','subject'})
% xticklabels({'shuffle-sound', 'sound', 'shuffle-category', 'category'});
ylabel('correct rate')
ylim([0 1])
title('Decoding with SVM')
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