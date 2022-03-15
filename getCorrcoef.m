function output = getCorrcoef(X, F, opt, varargin)

% usage:
% for component analysis
% opt_corr.CV = 1;
% opt_corr.comp = 1;
% output = getCorrcoef(X, F, opt_corr, Rs{1});

%% select features for regression
% In structure F:
% F_mat rows 1~9: frequency powers
% F_mat rows 10~10+7*9-1 = 10~72: combined spectrotemporal modulation
% powerd
% F_mat rows 73~end: full spectrotemporal modulation power (including
% negative and positive temporal rates)
% In structure C:
% for NatVoc - category_regressors: 1/0 labels for each sound
% for Nat (without voc) - continuous_Scores: 0~1, subject's rating

Y_freq = F.F_mat(1 : F.nFreq,:)';
Y_temp = F.F_mat(F.nFreq+1 : F.nFreq+F.nTemp,:)';
Y_spec = F.F_mat(F.nFreq+F.nTemp+1 : F.nFreq+F.nTemp+F.nSpec,:)';
Y_mod = F.F_mat(F.nFreq+F.nTemp+F.nSpec+1 : F.nFreq+F.nTemp+F.nSpec+F.nSpectemp,:)';
Y_mod_weighted = F.F_mat(F.nFreq+F.nTemp+F.nSpec+F.nSpectemp+1:end,:)';

% note: previous version of PCA with SVD did not center the data, didn't change the
% results much
% [U,S,V] = svd(Y_mod);
% Y_mod = U(:,1:15);
% [U,S,V] = svd(Y_mod_weighted);
% Y_mod_weighted = U(:,1:15); % 165*15
% figurex; plot(cumsum(diag(S)) ./ sum(diag(S)));

[~, Y_mod] = pca(Y_mod, 'NumComponents', 15);
[~, Y_mod_weighted] = pca(Y_mod_weighted, 'NumComponents', 15);

Y_cat = 100.*F.C.continuous_scores; % column number = predictor number, multiply by 100 to avoid regression error

% visualize PCs for spectromodulation
% figure
% for i = 1:15
%     V_pc = S(i,i).*V(:, i);
%     V_pc = reshape(V_pc, F.nSpec, F.nTemp);
%     subplot(3,5,i), imagesc(flipud(V_pc)), colormap(jet)
%     set(gca,'ytick',1:2:length(F.spec_mod_rates))
%     set(gca,'yticklabels',arrayfun(@num2str,fliplr(F.spec_mod_rates(1:2:end)),'UniformOutput',false))
%     set(gca,'xtick',1:2:length(F.temp_mod_rates)-1)
%     set(gca,'xticklabels',arrayfun(@num2str,F.temp_mod_rates(2:2:end),'UniformOutput',false))
%     axis square
%     set(gca,'fontsize',24);
% end

%% set parameters
% Y: features to be regress against (165 x N)
% z: response profile (165 x 1)

nStim = size(X,1);
nCase = 6; % possibilities of features combinations (freq, temp, spec, temp+spec, spectemp, cat)
N = size(F.F_mat,2); % number of sounds
fold = 5; % for 165 sounds: 11 fold --> 15 sounds; 5 fold --> 33 sounds

% initiate parameters
if opt.comp == 1 % for component analysis
    R = varargin{1};
    opt.K = size(R, 2);
else % for pixel-wise analysis
    opt.K = size(X,2);
    R = X;
end
zpredict = zeros(opt.K, N, nCase);
Coef = cell(1,nCase);
Rsquared = zeros(opt.K, nCase);
Corr = zeros(opt.K, nCase);

% plot order
% ind = [3 2 4 6 5 1]; % for 80Z
% ind = [1 2 5 4 3 6]; % for 132D 1/19 session
% ind = [1 2 4 6 5 3]; % for 132D, session 2
ind = 1:opt.K;

for i = 1:nCase
    switch i 
        case 1
            Y = Y_freq; 
        case 2
            Y = Y_temp; 
        case 3
            Y = Y_spec;      
        case 4
            Y = [Y_freq, Y_mod]; 
        case 5
            Y = [Y_freq, Y_temp, Y_spec]; 
        case 6
            Y = Y_cat; 
    end
    
    % === normalize features here ===
    Y = [ones(nStim,1), normalize(Y)];
%     Y = normalize(Y);
%      Y = [ones(nStim,1), Y];
    
    if ~opt.CV % no cross validation
        for k = 1:opt.K
            z = R(:, ind(k));
%             [b,~,~,~,stats] = regress(z, Y); 
%             % stats: R-squared (uncorrected), F-stat, p-value, error variance 
%             Rsquared(k,i) = stats(1);
%             Coef{i}(k,:) = b';
            output_regress = myRegression(z,Y,opt.regress_type);
            Rsquared(k,i) = output_regress.Rsquared;
            Coef{i}(k,:) = output_regress.coef;
           
        end
    else % cross validation
        seed = 42;
        rng(seed);
        c = cvpartition(N,'KFold',fold); % 11 fold, 15 sounds in each test group
        for k = 1:opt.K % components
            z = R(:, ind(k)); % responses
            for n = 1:fold % 11 folds
                cv_ind.test = find(test(c, n));
                cv_ind.training = find(training(c, n));                
                Ytest = Y(cv_ind.test, :);
                ztest = z(cv_ind.test);
                Ytrain = Y(cv_ind.training, :);     
                ztrain = z(cv_ind.training);
%                 [b,~,~,~,stats] = regress(ztrain, Ytrain);
                output_regress = myRegression(z,Y,opt.regress_type);
                b = output_regress.coef';
                zpredict(k, cv_ind.test, i) = Ytest*b;
                
                % measure 1.1: cv-Rsquared on each fold
                SSresid = sum((Ytest*b-ztest).^2);
                SStotal = (length(ztest) - 1) * var(ztest);
                rsq(n) = 1 - SSresid/SStotal; % uncorrected r-squared
            end
            cvR2_mean(k,i) = mean(rsq);
            cvR2_se(k,i) = std(rsq)/sqrt(length(rsq));
            
            % measure 1.2: cv-Rsquared on all datasets
%             SSresid = sum((zpredict(k, :, i)'-z).^2);
%             SStotal = (length(z) - 1) * var(z);
%             rsq = 1 - SSresid/SStotal; % uncorrected r-squared
            
            % measure 2: correlation
            rr = corrcoef(zpredict(k,:,i), z);
            Corr(k,i) = rr(1,2);

            % measure 3: normalized square error
%             Corr(k,i) = 1-getNSE(zpredict(k, :, i)', z, 0);

        end
        disp(['Case: ', num2str(i)])

    end
end

if ~opt.CV % if without cross validation
    output.Rsquared = Rsquared; 
    output.Coef = Coef;
else % if do cross validation
    output.zpredict = zpredict;
    output.Corr = Corr;
    output.cvR2_mean = cvR2_mean;
    output.cvR2_se = cvR2_se;
end

% plot figure
load('cmap_barplot_getCorrcoef_6colors.mat');
% cmap = zeros(nCase, 3);
% cmap(1,:) = [155, 192, 208]./255;
% cmap(2,:) = [3, 150, 166]./255;
% cmap(3,:) = [2, 83, 115]./255;
% cmap(4,:) = [242, 149, 68]./255;
% cmap(5,:) = [238, 213, 183]./255;
% cmap(6,:) = [150, 150, 150]./255;
    
if opt.comp
    if opt.CV
        % plot the correlation / NSE 
        figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w','position',[817 657 1224 464])
        hbar1 = bar(Corr);
        ylim([min(0, min(Corr(:))), 1])
        ylabel('Corr. Coef.')
%         ylabel('Normalized square error')
        legend({'Frequency',...
        'Temp Mod',...
        'Spec Mod',...
        'Frequency + SpecTemp Mod (joint)',...
        'Frequency + SpecTemp Mod (indep.)',...
        'Category'},'location','northeastoutside')
        xlabel('Component number')
        
        % plot cross-validated R^2
        figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w','position',[817 657 1224 464])
        hbar2 = bar(cvR2_mean); hold on
        er = errorbar(repelem(1:opt.K, size(cvR2_mean,2)) + repmat([hbar2.XOffset], 1, opt.K),... % positions or all bars
            reshape(cvR2_mean', numel(cvR2_mean),1),...
            reshape(cvR2_se', numel(cvR2_se),1));    
        er.Color = [0 0 0];                            
        er.LineStyle = 'none';  
%         ylim([min(0, min(cvR2_mean(:)-cvR2_se(:))) max(cvR2_mean(:)+cvR2_se(:))])
        ylim([0 1])
        ylabel('Cross-validated R^2')
        for ii = 1:nCase
            hbar1(ii).FaceColor = cmap(ii,:);
            hbar2(ii).FaceColor = cmap(ii,:);
        end
    else
        figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w','position',[817 657 1224 464])        
        hbar2 = bar(Rsquared);
        ylim([0 1])
        ylabel('Variance Explained (R^2)')
        for ii = 1:nCase
            hbar2(ii).FaceColor = cmap(ii,:);
        end
    end
    legend({'Frequency',...
        'Temp Mod',...
        'Spec Mod',...
        'Frequency + SpecTemp Mod (joint)',...
        'Frequency + SpecTemp Mod (indep.)',...
        'Category'},'location','northeastoutside')
    xlabel('Component number')
    
%     if ~opt.CV
%         % plot coefficients
%         figure, 
%         subplot(1,4,1), plot(Coef{1}')
%         subplot(1,4,2), plot(Coef{2}')
%         subplot(1,4,3), plot(Coef{3}')
%         subplot(1,4,4), imagesc(reshape(Coef{4}(3,10:end), F.nSpec, F.nTemp))
%     end
end

end
