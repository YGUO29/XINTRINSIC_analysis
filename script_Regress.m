% variance explained by acoustic features and categories
% Y: features to be regress against (165 x N)
% z: response profile (165 x 1)
% load('D:\SynologyDrive\=data=\F_halfcosine_marm_4reps.mat') 
load('D:\SynologyDrive\=data=\F_halfcosine_marm_NatVoc.mat') 

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
[U,S,V] = svd(Y_mod);
Y_mod = U(:,1:15);
[U,S,V] = svd(Y_mod_weighted);
Y_mod_weighted = U(:,1:15);

% visualize PCs for spectromodulation
% figure
% for i = 1:15
% V_pc = V(:, i);
% V_pc = reshape(V_pc, F.nSpec, F.nTemp);
% subplot(3,5,i), imagesc(flipud(V_pc)), colormap(jet)
% set(gca,'ytick',1:2:length(F.spec_mod_rates))
% set(gca,'yticklabels',arrayfun(@num2str,fliplr(F.spec_mod_rates(1:2:end)),'UniformOutput',false))
% set(gca,'xtick',1:2:length(F.temp_mod_rates)-1)
% set(gca,'xticklabels',arrayfun(@num2str,F.temp_mod_rates(2:2:end),'UniformOutput',false))
% axis square
% set(gca,'fontsize',24);
% end

Y_cat = 100.*F.C.continuous_scores; % column number = predictor number, multiply by 100 to avoid regression error

cmap = colormap(lines(7));

%% Variance explained - bar plot
% acoustics regression with components' response profiles
K = 6;
R = Rs{1};

nStim = size(X,1);
CV =0; % use cross validation or not
N = size(F.F_mat,2);
ztest = zeros(K, N);
zpredict = ztest;
Rsquared = ztest;
Corr = zeros(K,4);

% ind = [3 2 4 6 5 1]; % for 80Z
% ind = [1 2 5 4 3 6]; % for 132D 1/19 session
% ind = [1 2 4 6 5 3]; % for 132D, session 2
ind = 1:K;

for i = 1:4
    switch i 
        case 1
            Y = [ones(nStim,1), Y_freq]; 
        case 2
            Y = [ones(nStim,1), Y_temp]; 
        case 3
            Y = [ones(nStim,1), Y_spec];      
        case 4
%             Y = [ones(nStim,1), Y_freq, Y_mod];           
            Y = [ones(nStim,1), Y_freq, Y_temp, Y_spec]; 
    end
    
    if ~CV
        for k = 1:K
            z = R(:, ind(k));
            [b,~,~,~,stats] = regress(z, Y); 
            % stats: R-squared, F-stat, p-value, error variance
            Corr(k,i) = stats(1);
        end
    else
        % c = cvpartition(N,'HoldOut',1);
        for k = 1:K
            z = R(:, ind(k));
            for n = 1:N % leave-one out, use the rest for training 
                c.test = n;
                c.training = setdiff(1:N, n);
            % ztest = zeros(K, 1);
            % zpredict = ztest;
            % Rsquared = ztest;
                Ytest = Y(c.test, :);
                Ytrain = Y(c.training, :);
        %         ztest(k) = z(c.test);
                ztrain = z(c.training);
                [b,~,~,~,stats] = regress(ztrain, Ytrain); 
                Rsquared(k, n) = stats(1);
                zpredict(k, n) = Ytest*b;
            end
            rr = corrcoef(zpredict(k,:), z);
            Corr(k,i) = rr(1,2).^2;
        end
    end
end
figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w','position',[817   701   857   420])
hbar = bar(Corr)
hbar(1).FaceColor = cmap(4,:);
hbar(2).FaceColor = cmap(5,:);
hbar(3).FaceColor = cmap(6,:);
hbar(4).FaceColor = cmap(2,:);

legend({'Frequency only',...
    'Temp Mod',...
    'Spec Mod',...
    'All acoustic'},'location','northeastoutside')
ylabel('Variance Explained')
xlabel('Component number')

% xticklabels({'Freq',...
%     'Freq+Temp',...
%     'Freq+Spec',...
%     'All'})
% xtickangle(45)

%% Variance explained - bar plot
% acoustics & category regression with components' response profiles

CV = 0; %use cross validation or not
N = size(F.F_mat,2);
ztest = zeros(K, N);
zpredict = ztest;
Rsquared = ztest;
Corr = zeros(K,3);

% ind = [3 2 4 6 5 1]; % for 80Z
% ind = [1 2 5 4 3 6]; % for 132D 1/19 session
% ind = [1 2 4 6 5 3]; % for 132D, session 2
ind = 1:K;
% Rsquared = zeros(K, 3);
for i = 1:3
    switch i 
        case 1
            Y = [ones(nStim,1), Y_cat]; 
        case 2
            Y = [ones(nStim,1), Y_freq, Y_mod]; 
        case 3
            Y = [ones(nStim,1), Y_freq, Y_mod, Y_cat];           
    end
 
    if ~CV
        for k = 1:K
        z = R(:, ind(k));
        [b,~,~,~,stats] = regress(z, Y); 
        % stats: R-squared, F-stat, p-value, error variance
        Corr(k,i) = stats(1);
        end
    else
        for k = 1:K
            z = R(:, ind(k));
            for n = 1:N % leave-one out, use the rest for training 
                c.test = n;
                c.training = setdiff(1:N, n);
            % ztest = zeros(K, 1);
            % zpredict = ztest;
            % Rsquared = ztest;
                Ytest = Y(c.test, :);
                Ytrain = Y(c.training, :);
        %         ztest(k) = z(c.test);
                ztrain = z(c.training);
                [b,~,~,~,stats] = regress(ztrain, Ytrain); 
                Rsquared(k, n) = stats(1);
                zpredict(k, n) = Ytest*b;
            end
            rr = corrcoef(zpredict(k,:), z);
            Corr(k,i) = rr(1,2);
        end
    end
end
figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w','position',[817   701   857   420])
hbar = bar(Corr)
legend({'Category only',...
    'All Acoustic',...
    'Acoustic + Category'},'location','northeastoutside')
ylabel('Variance Explained')
xlabel('Component number')

% xticklabels({'Category',...
%     'All Acoustics',...
%     'Category+Acoustics'})
% xtickangle(45)

%% Calculate R-square 
% regression with acoustic & category for each pixel
CV = 0;
nPix = size(X,2);
nStim = size(X,1);
ind = 1:nPix;

Corr = [];
Coef = cell(1,6);
for i = 1:6
    switch i 
        case 1
%             Y = [ones(size(Y_cat)), Y_temp, Y_spec]; % delete feature of interest
            Y = [ones(nStim,1), Y_freq]; 
        case 2
%             Y = [ones(size(Y_cat)), Y_freq, Y_spec]; 
            Y = [ones(nStim,1), Y_temp]; 
        case 3
%             Y = [ones(size(Y_cat)), Y_freq, Y_temp];      
            Y = [ones(nStim,1), Y_spec];      
        case 4
            Y = [ones(nStim,1), Y_freq, Y_mod];  
        case 5
            Y = [ones(nStim,1), Y_freq, Y_temp, Y_spec];
        case 6
            Y = [ones(nStim,1), Y_cat]; 
    end
    if CV
    % c = cvpartition(N,'HoldOut',1);
        for k = 1:nPix
            z = X(:, ind(k));
            for n = 1:nStim % leave-one out, use the rest for training 
                c.test = n;
                c.training = setdiff(1:nStim, n);
            % ztest = zeros(K, 1);
            % zpredict = ztest;
            % Rsquared = ztest;
                Ytest = Y(c.test, :);
                Ytrain = Y(c.training, :);
        %         ztest(k) = z(c.test);
                ztrain = z(c.training);
                [b,~,~,~,stats] = regress(ztrain, Ytrain); 
                Rsquared(k, n) = stats(1);
                zpredict(k, n) = Ytest*b;
            end
            rr = corrcoef(zpredict(k,:), z);
            Corr(k,i) = rr(1,2);
        end
    else
        for k = 1:nPix
            z = X(:, ind(k));
            
            % ==== regular regression ====
%             [b,~,~,~,stats] = regress(z, Y); 
%             % stats: R-squared, F-stat, p-value, error variance
%             Corr(k,i) = stats(1);
%             Coef{i}(k,:) = b';

            % ==== lasso regression ====
            [B, FitInfo] = lasso(Y, z,'CV', 3);
            b = B(:, FitInfo.IndexMinMSE); % IndexMinMSE or Index1SE
            b0 = FitInfo.Intercept(FitInfo.IndexMinMSE);
            zpredict = Y*b + b0;
            rr = corrcoef(zpredict, z);
            Corr(k,i) = rr(1,2);
            Coef{i}(k,:) = b';
%             lassoPlot(B,FitInfo,'PlotType','CV');
%             legend('show') % Show legend
      
            % for checking residuals
%             if k == 4500
%             z_hat = Y*b;
%             figurex; scatter(z_hat, (z_hat - z));
%             end
        end
    end
i
end

% make difference maps
% for i = [1:3, 5]
%     Corr(:,i) = Corr(:,4) - Corr(:,i);
% end

%% coefficient map

Title = {'Frequency (Hz)',...
    'Temp. Mod. (Hz)',...
    'Spec. Mod. (Cyc/Oct)',...
    'Frequency and SpecTemp. Mod.',...
    'Category',...
    'Category (adjusted colorbar)'};
% Title = {'Frequency',...
%     'Frequency and Temp. Mod.',...
%     'Frequency and Spec. Mod.',...
%     'Frequency and SpecTemp. Mod.',...
%     'Category',...
%     'Category (adjusted colorbar)'};
figurex([563         253        1565         883])

for i = [1:3, 6]
    % calculated weighted average of features (the feature with highest
    % coefficients)
    b = Coef{i}(:,2:end); % 9000 x nFeature
    corr = Corr(:,i);
    switch i
        case 1 % frequency
            feat = F.FreqBounds(2:end); % upper bounds
            cb_labels = arrayfun(@num2str,feat,'UniformOutput',false);
        case 2 % spectral modulation
            feat = F.temp_mod_rates(2:end); 
            cb_labels = arrayfun(@num2str,feat,'UniformOutput',false);
        case 3 % temporal modulation
            feat = F.spec_mod_rates;
            cb_labels = arrayfun(@num2str,feat,'UniformOutput',false);
        case 6 % category
            cb_labels = F.C.category_labels;
            feat = 1:length(cb_labels);

    end
      
    % ==== weighted average as a visualization measurement (negative coefficients problem!!)
%     b = -b;
%     b(b<0) = 0;
%     b = b./repmat(sum(b,2), 1, size(b,2)); 
%     temp = b*[1:length(feat)]'; 
    % ==== max as a visualization measurement ====
    [~,b_ind] = max(b,[],2); 

    b_ind = reshape(b_ind,[para.height, para.width]);
    map = repmat(b_ind, 1, 1, 3);
    map(:,:,1) = (map(:,:,1)./length(feat)).*0.8; % H
%     map(:,:,1) = (map(:,:,1)./max(max(map(:,:,1)))).*0.8; % H
    map(:,:,2) = reshape(corr,[para.height, para.width]); % S
    map(:,:,2) = map(:,:,2)./max(max(map(:,:,2))); % S
    map(:,:,3) = map(:,:,2); % V
%     map(:,:,2) = ones(size(map(:,:,2)));
    map = hsv2rgb(map);
    
    if i == 6
    subplot(2,2,4)
    else
    subplot(2,2,i)
    end
    imagesc(map), axis image, 
    
    cmap_hsv = ones(length(feat), 3);
    cmap_hsv(:,1) = ((1:length(feat))./length(feat)).*0.8';
    cmap = hsv2rgb(cmap_hsv);
%     cmap = hsv(length(feat));
    
    colormap(gca, cmap);
    hcb(i) = colorbar;
    cb_ticks = range(hcb(i).Limits)/length(feat):range(hcb(i).Limits)/length(feat):hcb(i).Limits(2);
    set(hcb(i), 'Ticks', cb_ticks, 'TickLabels',cb_labels);
    title(Title{i})
    axis off
    
%     plotContour(ct);

end
%% R-squared map
Corr1 = Corr(:,1:6);
Corr2 = Corr(:,6);
Title = {'Frequency',...
    'Temp. Mod.',...
    'Spec. Mod.',...
    'Frequency + SpecTemp. Mod.',...
    'Frequency + Spec + Temp. Mod.',...
    'Category',...
    'Category (adjusted colorbar)'};
% Title = {'Frequency',...
%     'Frequency and Temp. Mod.',...
%     'Frequency and Spec. Mod.',...
%     'Frequency and SpecTemp. Mod.',...
%     'Category',...
%     'Category (adjusted colorbar)'};
figurex([563         253        1565         553])
Max = max(Corr1(:));
% Max = 0.3;
Min = prctile(Corr1(:),0);
for i = 1:6
    subplot(2,3,i)
    b_ind = reshape(Corr1(:,i),[para.height, para.width]);
%     Min = prctile(temp(:),75);
    Min = 0;
    imagesc(b_ind, [Min, Max]), axis image, colorbar
    title(Title{i})
    axis off
    colormap(jet)
%     plotContour(ct);

end

figurex([1031         711        1756         627])
Max = max(Corr(:));
Min = prctile(Corr2(:),0);
for i = 1:2
    subplot(1,2,i)
    if i== 2
        b_ind = reshape(Corr(:,6),[para.height, para.width]);
        Max = max(b_ind(:));
%         Min = prctile(Corr(:),0);
        Min = 0;
    else
        b_ind = reshape(Corr(:,6),[para.height, para.width]);
    end
%     Max = max(temp(:));
    imagesc(b_ind, [Min, Max]), axis image, colorbar
    title(Title{i+5})
    axis off
    colormap(jet)
%     plotContour(ct);

end

%% overlay normalized regression map for frequency, temp mod, spec mod in RGB channel
figurex([669,738,1263,420]); 
img_rgb = zeros(para.height, para.width, 3);
for i = 1:3
    b_ind = reshape(Corr(:,i),[para.height, para.width]);
    b_ind = b_ind-min(b_ind(:));
%     temp = temp./max(temp(:));
    img_rgb(:,:,i) = b_ind;
end
img_rgb = img_rgb./max(img_rgb(:));
subplot(1,2,1)
imagesc(img_rgb), axis image, axis off
title('R=freq, G=temp.mod., B=spec.mod.')

% overlay MAX of regression map for frequency, temp mod, spec mod in RGB channel
img_rgb = zeros(para.height, para.width, 3);
b_ind = reshape(Corr(:,1:3),[para.height, para.width, 3]);
[~, ind_max] = max(b_ind,[],3);
for i = 1:para.height
    for j = 1:para.width
%         img_rgb(i,j,ind_max(i,j)) = 1;
        img_rgb(i,j,ind_max(i,j)) = b_ind(i,j, ind_max(i,j));
    end
end
img_rgb = (img_rgb - min(img_rgb(:)))./range(img_rgb(:));
subplot(1,2,2)
imagesc(img_rgb), axis image, axis off
title('R=freq, G=temp.mod., B=spec.mod.')