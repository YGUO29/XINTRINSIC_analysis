% perform DSS
%% whiten single rep data matrix
% DataMat = DataMat(:,:,6:95,:,:); para.height = 90;
opt = struct;
[X, mov_rel, mov_rel_sep] = getX(DataMat, para, opt);

%% divide into train and test sets
% averaged response (dur stimulus presentation time window)
D = squeeze(mean(mov_rel_sep(:,:,:,:,para.fr*para.preStim+1:para.fr*(para.preStim+para.durStim)),5));

% ==== for natural + model-matched sounds ====
% D_ori = D(:,1:2:para.nStim-1,:,:); % original
% D_ori = reshape(D_ori, para.nRep, para.nStim/2, para.width*para.height);
% D_mm = D(:,2:2:para.nStim,:,:); % model-matched
% D_mm = reshape(D_mm, para.nRep,  para.nStim/2, para.width*para.height);
% D_all_3d = reshape(D, para.nRep, para.nStim, para.width*para.height);
% 
% seed_train = randperm(para.nStim/2, para.nStim*3/8);
% seed_test = setdiff(1:para.nStim/2, seed_train);
% D_train = cat(2, D_ori(:,seed_train,:), D_mm(:,seed_train,:));
% D_test = cat(2, D_ori(:,seed_test,:), D_mm(:,seed_test,:));

% ==== for natural sounds only ====
% D_all_3d = reshape(D, para.nRep, para.nStim, para.width*para.height);
% seed_train = randperm(para.nStim, para.nStim*4/5);
% seed_test = setdiff(1:para.nStim, seed_train);
% D_train = D_all_3d(:,seed_train,:);
% D_test = D_all_3d(:,seed_test,:);

% ==== for vocalizations ====
% D_ori = D(:,1:2:para.nStim-1,:,:); % original
% D_ori = reshape(D_ori, para.nRep, para.nStim/2, para.width*para.height);
% D_mm = D(:,2:2:para.nStim,:,:); % model-matched
% D_mm = reshape(D_mm, para.nRep,  para.nStim/2, para.width*para.height);
D_all_3d = reshape(D, para.nRep, para.nStim, para.width*para.height);
% 
seed_train = randperm(para.nStim/4, para.nStim*3/16);
seed_test = setdiff(1:para.nStim/4, seed_train);
D_train = D_all_3d(:,[seed_train, seed_train+para.nStim/4, seed_train+para.nStim/2, seed_train+para.nStim*3/4], :);
D_test = D_all_3d(:,[seed_test, seed_test+para.nStim/4, seed_test+para.nStim/2, seed_test+para.nStim*3/4], :);

% D_avg       = cat(1, squeeze(mean(D_ori,1)), squeeze(mean(D_mm,1)));
D_avg       = squeeze(mean(D_all_3d, 1));
D_train_avg = squeeze(mean(D_train, 1));
D_test_rep1 = squeeze(mean(D_test(1:2:para.nRep, :,:), 1));
D_test_rep2 = squeeze(mean(D_test(2:2:para.nRep, :,:), 1));
% D_test_rep1 = squeeze(mean(D_test(1:para.nRep/2, :,:), 1));
% D_test_rep2 = squeeze(mean(D_test(para.nRep/2+1:end, :,:), 1));
D_test_allrep = squeeze(mean(D_test,1,'omitnan'));

% de-mean data
for i = 1:size(D_train_avg,2)
    D_avg(:,i) = D_avg(:,i) - mean(D_avg(:,i));
    D_train_avg(:,i) = D_train_avg(:,i) - mean(D_train_avg(:,i));
    D_test_rep1(:,i) = D_test_rep1(:,i) - mean(D_test_rep1(:,i));
    D_test_rep2(:,i) = D_test_rep2(:,i) - mean(D_test_rep2(:,i));
    D_test_allrep(:,i) = D_test_allrep(:,i) - mean(D_test_allrep(:,i));
end
%% perform DSS (for cross-validation)
 
% whiten data matrix
D_whiten = zeros(size(D_train));
for iRep = 1:para.nRep
    D_temp = squeeze(D_train(iRep,:,:));
    [U, S, V] = svd(D_temp, 'econ');
    D_whiten(iRep, :, :) = U*V';
%     D_whiten(iRep, :, :) = (U(:,1:size(V,1))*V')';
end

% average whitened matrix across repetitions
D_whiten_avg = squeeze(mean(D_whiten,1));

%% cross-validation to pick number of components
nK = 20;
corr_matrix = zeros(para.width*para.height, nK);
% PCA on averaged data matrix
[U, S, V] = svd(D_whiten_avg,'econ');

for k = 1:nK
R = U(:,1:k)*S(1:k,1:k);
% W = V(:,1:k)';
W = pinv(R)*D_train_avg;

% % project onto test dataset
R_test = D_test_rep1 * pinv(W);
D_test_rep1_denoised = R_test * W;
R_test = D_test_rep2 * pinv(W);
D_test_rep2_denoised = R_test * W;
% 
% % calculate correlation
% corr_matrix(:,k) = diag(corr(D_test_rep1_denoised, D_test_rep2));

R_test_allrep = D_test_allrep * pinv(W);
D_test_allrep_denoised = R_test_allrep * W;
corr_matrix(:,k) = diag(corr(D_test_allrep_denoised, D_test_allrep));
end

% plot upper-bound and one-split correlation of each pixel
upper_bound = diag(corr(D_test_rep1_denoised, D_test_rep2_denoised));
figurex;
subplot(1,2,1)
imagesc(reshape(upper_bound, para.height, para.width),[0 1]), axis image, title('upper bound')
colorbar
subplot(1,2,2)
imagesc(reshape(corr_matrix(:,6), para.height, para.width), [0 1]), axis image, title('split-corr')
colorbar

figurex;
SEM     = std(corr_matrix)/sqrt(size(corr_matrix,1));               % Standard Error
ts      = tinv([0.025  0.975],size(corr_matrix,1)-1);      % T-Score
shadedErrorBar(1:nK, median(corr_matrix, 'omitnan'),ts(2).*SEM,'lineProps','b');
hold on, plot([6 6], [0.2 1],'--')
% upper bound
% hold on, 
% SEM     = std(upper_bound)/sqrt(length(upper_bound));               % Standard Error
% ts      = tinv([0.025  0.975],length(upper_bound) -1);      % T-Score
% shadedErrorBar(1:nK, ones(1,nK).*median(upper_bound),ones(1,nK).*ts(2).*SEM,'lineProps','c');
xlabel('#components')
ylabel('correlation before & after component-projection')

%% plot scatter plots
% for i = 1:size(D_train_avg,2)
%     D_test_rep1(:,i) = D_test_rep1(:,i) - mean(D_test_rep1(:,i));
%     D_test_rep2(:,i) = D_test_rep2(:,i) - mean(D_test_rep2(:,i));
%     D_test_rep1_denoised(:,i) = D_test_rep1_denoised(:,i) - mean(D_test_rep1_denoised(:,i));
%     D_test_rep2_denoised(:,i) = D_test_rep2_denoised(:,i) - mean(D_test_rep2_denoised(:,i));
% end
figurex;

% hist3([x,y],'Nbins',[100 100],'CdataMode','auto'), view(2)
x = diag(corr(D_test_rep1, D_test_rep2));
y = diag(corr(D_test_rep2_denoised, D_test_rep1_denoised));
hold on, scatter(x,y,'.')

x = diag(corr(D_test_rep1, D_test_rep2));
y = diag(corr(D_test_rep2, D_test_rep1_denoised));
hold on, scatter(x,y,'.')

axis([-0.2 1 -0.2 1])
hold on, plot([-0.2, 1], [-0.2, 1]);
hold on, plot(0:0.01:0.8, sqrt(0:0.01:0.8));
%% denoise using R
% whiten the entire dataset
D_whiten = zeros(size(D_all_3d));
for iRep = 1:para.nRep
    D_temp = squeeze(D_all_3d(iRep,:,:));
    [U, S, V] = svd(D_temp, 'econ');
    D_whiten(iRep, :, :) = U*V';
%     D_whiten(iRep, :, :) = (U(:,1:size(V,1))*V')';
end
% average whitened matrix across repetitions
D_whiten_avg = squeeze(mean(D_whiten,1));

k = 8;
[U, S, V] = svd(D_whiten_avg);
R = U(:,1:k)*S(1:k,1:k);

X_denoise = R*pinv(R)*X;

% ==== visualize denoised images vs. original images ====
figurex;
for i = 1:144
    subplot(1,2,1), imagesc(reshape(X(i,:), para.height, para.width)), axis image;
    title(['Original, ', C.category_labels(C.category_assignments(i))])
    subplot(1,2,2), imagesc(reshape(X_denoise(i,:), para.height, para.width)), axis image;
    title(['Denoised, ', C.category_labels(C.category_assignments(i))])
    pause
end

% run decomposition in MATLAB
opt.fluo = 1; 
opt.method = 'mICA'; % 'mICA' or 'NMF' or 'PCA'
opt.nRows = 1;
opt.plotON = 1;
Ks = k;
% opt.X_test = X_test;
[Rs, Ws, comp, recon_error, X_hats] = runDecomp(X_denoise, Ks, opt, para);