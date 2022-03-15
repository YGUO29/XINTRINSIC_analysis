% Compute variance explained by acoustic features and categories
% load('D:\SynologyDrive\=data=\F_halfcosine_marm_4reps.mat') 

% clear X
% load('E:\Dropbox\_Yueqi Guo\_research\Imaging\=data=\102D\DataMatProc_Calcium_210226_102D_NatVocMM_12reps.mat')
ind_voc_orig = 1:2:30; ind_voc_synth = 2:2:30; 
ind_nat_orig = 31:2:size(X,1); ind_nat_synth = 32:2:size(X,1); 
X_nat_orig = X(ind_nat_orig, :);
X_voc_orig = X(ind_voc_orig, :);
X_voc_synth = X(ind_voc_synth, :);
X_nat_synth = X(ind_nat_synth, :);

%% Component analysis: Variance explained for each component - bar plot
% acoustics regression with components' response profiles
% ==== run decomposition first ====
opt.fluo = 1; 
opt.method = 'mICA'; % 'mICA' or 'NMF' or 'PCA'
opt.nRows = 1;
opt.plotON = 1;
opt.reverse = 1;
Ks = 6;
% opt.X_test = X_test;
figurex;
[Rs, Ws, comp, recon_error, X_hats] = runDecomp(X, Ks, opt, para);

load('D:\SynologyDrive\=data=\F_halfcosine_marm_NatJM.mat') 
% load('D:\SynologyDrive\=data=\F_halfcosine_marm_NatVoc.mat')
opt_corr.CV = 1;
opt_corr.comp = 1;
opt_corr.regress_type = 'regular';
output = getCorrcoef(X, F, opt_corr, Decomp_132D.Rs{2});

%% predict response to other sounds

for iCase = 1:6

% X_diff = abs(mean((X_voc_orig - X_hat)./X_voc_orig, 1));
X_diff = abs(mean(X_true - X_hat, 1));

figurex; 
subplot(1,3,1)
imagesc(reshape(mean(X_true,1), para.height, para.width),[0 0.22]), axis image, axis off
colorbar
subplot(1,3,2)
imagesc(reshape(mean(X_hat,1), para.height, para.width),[0 0.22]), axis image, axis off
colorbar
subplot(1,3,3)
imagesc(reshape(X_diff, para.height, para.width)), 
axis image, axis off
%     colormap(cbrewer('div', 'PRGn',256))
colormap(hot)
colorbar

end

%% reconstruct X with predicted Rs
clear X_hat
% load('D:\SynologyDrive\=data=\F_halfcosine_marm_NatVoc.mat')
% predict vocalizations
X_true = X_voc_synth;
feature = F.table(1:15, :);
nStim = size(X_true, 1);
nComp = size(Ws{1}, 1);
% predict synthesized sounds
% X_true = X_nat_synth;
% feature = F.table(16:end, :);
nCase = 5;

Titles = {'Frequency',...
        'Temp Mod',...
        'Spec Mod',...
        'Frequency + SpecTemp Mod (joint)',...
        'Frequency + SpecTemp Mod (indep.)',...
        'Category'};
    
Recon_error = zeros(nStim, nCase);
Corr_pix = Recon_error;
Corr_comp = zeros(nComp, nCase);
for iCase = 1:nCase
    
    % ==== estimate X_hat: from CV output ====
%     if opt_corr.comp % for component-based analysis:
%         X_hat = squeeze(output.zpredict(:,:,iCase))'*Ws{1};
%     else% for pixel-based analysis:
%         X_hat = squeeze(output.zpredict(:,:,iCase))';
%     end
    
    % ==== estimate X_hat: from coefficients ====

    switch iCase
        case 1
            feature_temp = feature.freq'; % #feature x #stim
        case 2
            feature_temp = feature.temp';
        case 3
            feature_temp = feature.spec';
        case 4 % project spectemp features to the pca space obtained from natural sounds
%             Y_mod = F.table.spectemp_mean(16:end, :);
            Y_mod = F.table.spectemp_mean;
            coeff = pca(Y_mod, 'NumComponents', 15); % 
            feature_temp = feature.spectemp_mean - mean(feature.spectemp_mean,1);
            feature_temp = [feature.freq'; (feature_temp * coeff)']; % #feature x #stim
        case 5
            feature_temp = [feature.freq';feature.temp';feature.spec'];
            
        otherwise
    end

    % calculate measurement for 'R' component responses
    R_hat = output.Coef{iCase}*[ones(1,size(feature_temp, 2));feature_temp]; % [#comp x #feat] x [#feat x #stim]
    R_hat = R_hat'; % #stim x #components
    R_proj = X_true * pinv(Ws{1});
    for iComp = 1:nComp
        rr = corrcoef(R_hat(:, iComp), R_proj(:, iComp));
        Corr_comp(iComp, iCase) = rr(1, 2);
%         SSresid = sum((R_hat(:, iComp) - R_proj(:, iComp)).^2);
%         SStotal = (length(R_proj(:, iComp)-1)) * var(R_proj(:, iComp));
%         Corr_comp(iComp, iCase) = 1 - SSresid/SStotal;
    end
    
    % calculate measurement for response patterns,
    X_hat(:,:,iCase) = R_hat*Ws{1};
    for iStim = 1:nStim
    %     subplot(1,2,1)
    %     imagesc(reshape(X(i,:), para.height, para.width)), axis image
    %     subplot(1,2,2)
    %     imagesc(reshape(X_hat(i,:), para.height, para.width)), axis image
    %     pause
        recon_error = norm(X_true(iStim,:) - X_hat(iStim,:,iCase), 'fro');
        corr_temp = corrcoef(X_true(iStim,:), X_hat(iStim,:,iCase));
        Recon_error(iStim, iCase) = recon_error; 
        Corr_pix(iStim, iCase) = corr_temp(1,2);
    end
   

    X_diff = abs(mean(X_true - X_hat(:,:,iCase), 1));
    figurex; 
    imagesc(reshape(X_diff, para.height, para.width)), 
    axis image, axis off
%     colormap(cbrewer('div', 'PRGn',256))
    colormap(hot)
    colorbar
    title(Titles{iCase});
end

% plot explained component response
figurex([1440         918        1145         420]);
hbar = bar(Corr_comp);
load('cmap_barplot_getCorrcoef_6colors.mat');
cmap = cmap(1:5, :);
for ii = 1:nCase
    hbar(ii).FaceColor = cmap(ii,:);
end
ylim([min(0, min(Corr_comp(:))) 1])
title('predicted & projected component responses')
legend({'Frequency',...
    'Temp Mod',...
    'Spec Mod',...
    'Frequency + SpecTemp Mod (joint)',...
    'Frequency + SpecTemp Mod (indep.)',...
    'Category'},'location','northeastoutside')
ylabel('Corr. Coef.')
xlabel('Component number')


% plot pixel-based correlation and reconstruction error
figure('color', 'w'); 
yyaxis left, plot(mean(Recon_error), '-*'), ylabel('recon error')
yyaxis right, plot(mean(Corr_pix), '-o'), ylabel('correlation')
xticks(1:6)
xticklabels(Titles)
xtickangle(45)
[mean(Recon_error); mean(Corr_pix)]
%% visualize predicted & actual response

figurex; 
for i = 1:nStim  
    subplot(2,4,1), imagesc(reshape(X_hat(i,:,1), ...
        para.height, para.width)), axis image
    title('predicted - frequency')
    subplot(2,4,2), imagesc(reshape(X_hat(i,:,2), ...
        para.height, para.width)), axis image
    title('predicted - temp.mod.')
    subplot(2,4,3), imagesc(reshape(X_hat(i,:,3), ...
        para.height, para.width)), axis image
    title('predicted - spec.mod.')
    subplot(2,4,4), imagesc(reshape(X_hat(i,:,4), ...
        para.height, para.width)), axis image
    title('predicted - spectemp.mod.')
    subplot(2,4,5), imagesc(reshape(X_hat(i,:,5), ...
        para.height, para.width)), axis image
    title('predicted - spec. + temp. mod.')
%     subplot(2,4,6), imagesc(reshape(X_hat(i,:,1), ...
%         para.height, para.width)), axis image
%     title('predicted - category')
    subplot(2,4,7), imagesc(reshape(X_true(i,:), ...
        para.height, para.width)), axis image; 
    title('actual response')
    pause; 
end

