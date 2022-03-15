function X_pred = getPrediction(X, para, K)
% this method is wrong...
%%
feature_file = uigetfile('D:\SynologyDrive\=data=\F*.mat');
load(fullfile('D:\SynologyDrive\=data=\', feature_file))

for i = 1:K
    for j = 1:165
        beta(i, j) = corr2(F.spectemp_mod(:,:,j), spectemp_r{i});
    end
end
X_pred = X;
for iSound = 1:165
    X_pred(iSound,:) = beta(:,iSound)'*Ws{1};
end
%%
figure,

for iSound = 1:165
subplot(1,2,1)
imagesc(reshape(X(iSound,:), para.height, para.width)), axis image
subplot(1,2,2)
imagesc(reshape(X_pred(iSound,:), para.height, para.width)), axis image
pause
end
end