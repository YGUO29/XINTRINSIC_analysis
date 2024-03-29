% Component analysis including reliability analysis

%% ============= Reliability test 1, select Perc of the trials/reps ==========
%% 
reduced_variable = 'trial'; % rep or trial
% Reps_temp = zeros(nRep, 4);
% for i = 1:nRep
%   Reps_temp(i,:) =  datasample(1:para.nRep, 4, 'replace',false);
% end

% change #selected trials, by randomly select certain percent
Perc = [fliplr(0.7:0.1:1), 0.5, 0.3, 0.1, 0.035];
Rep = [1:6]/6;
if strcmp(reduced_variable,'trial')
    rv = Perc;
elseif strcmp(reduced_variable, 'rep')
    rv = Rep;
else
end
nIter = 5;
% W2 = cell(length(rv), nIter);
% for i = length(Reps)
for i = 1
    for k = 1:nIter % repeat for nRep times
        if strcmp(reduced_variable,'trial')
            trials =    datasample(1:para.nStim,round(rv(i)*para.nStim),'replace',false); 
            trials = sort(trials,'ascend');
    %         opt.reps        = Reps_temp(k,:);
    %         opt.reps = datasample(1:para.nRep, Reps(i), 'replace',false)
            X_temp = squeeze( mean(X_sep(:, trials, :), 1) );        
            disp(['ReducedVariable=',num2str(i),', Iter=',num2str(k),' Selected ',num2str(round(rv(i)*para.nStim)),' trials'])

        elseif strcmp(reduced_variable, 'rep')
            reps =    datasample(1:para.nRep, round(rv(i)*para.nRep), 'replace',false); 
%         opt.reps        = Reps_temp(k,:);
%         opt.reps = datasample(1:para.nRep, Reps(i), 'replace',false)
            X_temp = squeeze( mean(X_sep(reps, :, :), 1) );
            disp(['ReducedVariable=',num2str(i),', Iter=',num2str(k),' Selected ',num2str(round(rv(i)*para.nRep)),' repetitions'])

        else
            disp('set reduced_variable as rep or trial')
        end
        % ======== perform ICA analysis =============
        opt_decomp.fluo = 1; 
        opt_decomp.method = 'mICA'; % 'mICA' or 'NMF' or 'PCA'
        opt_decomp.nRows = 1;
        opt_decomp.plotON = 0;
        opt_decomp.Ks = 6;
        para.mirror = 1;
        figurex;
        Decomp_temp = runDecomp(X_temp, opt_decomp, para);
        W2{i,k} = Decomp_temp.Ws{1};

    end
end


%% find match with original weights
result.dist = zeros(length(rv), nIter); % minimal distance
result.ind = cell(length(rv), nIter);
P = perms(1:opt_decomp.Ks); % all possible orders
% figure
% set(gcf,'color','w','position', [1323, 380, 1365, 192]);
for m = 1:length(rv)
    for n = 1:nIter
        Dist = zeros(1,size(P,1)); % length = number of possible permutations
        W_temp = W2{m,n};
        for j = 1:size(P,1) % go over all possible permutations for W2
            ind = P(j,:);
            W_temp2 = W_temp(ind,:); % rearranged W_temp
%             Dist(j) = sum(vecnorm(W_temp2' - Decomp_102D.Ws{1}'));
            Dist(j) = sum(vecnorm(abs(W_temp2)' - abs(Decomp_102D.Ws{1})'));
        end
        [result.dist(m,n), indDist] = min(Dist);
        result.ind{m,n} = P(indDist,:);

        % W1 = W(P(indDist,:),:);

        % ======== plot ============
        
%         set(gcf,'color','w','position', [1323, 380, 1365, 192]);
%         p = [1 6];
%         ha = tight_subplot(p(1),p(2),[0 0],[.1 .01],[.01 .01]);
%         for i = 1:K
%             cutoff = 0.2;
%             comp{i} = zeros(para.height, para.width);
%             comp{i}(ind_save) = W_temp(ind(i),:);
%             axes(ha(i)); imagesc(comp{i},cutoff.*[-1 1]),axis image, colormap(jet)
%             axis off
%         end
%         pause
    end
end
%% plot the best matched 'distance'
yy = mean(result.dist,2);
ee = std(result.dist,[],2)./sqrt(size(result.dist,2));
figurex([1440         802         560         536]);
if strcmp(reduced_variable,'trial')
    xx =rv.*para.nStim;
    errorbar(xx, yy, ee,'linewidth',2, 'Color', 'k')
    xlabel('Number of trials selected')
    ylabel('Mean Euclidean distance of components')
    xlim([0 para.nStim])
elseif strcmp(reduced_variable,'rep')
    xx = rv.*para.nRep;
    errorbar(xx, yy, ee,'linewidth',2, 'Color', 'k')
    xlabel('Number of repetitions selected')
    ylabel('Mean Euclidean distance of components')
    xlim([2 para.nRep-2])
end
hold on, plot([min(xx), max(xx)], ...
    (min(yy)+max(yy))./2.*ones(1,2), ...
    'color', 0.5.*ones(1,3), ...
    'linestyle', '--')
axis square

interp1(yy, xx, (min(yy)+max(yy))./2)


%% plot new components with a certain Perc number

for i = 1:nIter
W_temp = W2{8, i};
W_temp = W_temp(result.ind{8,i},:);

comp_temp = [];
figurex([1073         539         560          64]);
ha = tight_subplot(1, size(W_temp,1), [.01 .01],[.1 .01],[.01 .01]);
    for j = 1:size(W_temp,1)
        axes(ha(j)); 
        comp_temp = reshape(W_temp(j,:), para.height, para.width);
        cutoff = mean(W_temp(j,:)) + 7*std(W_temp(j,:));
        imagesc(comp_temp,  cutoff.*[-1 1]), axis image, colormap(jet), axis off
        drawnow
        if para.mirror
            set(gca,'XDir','reverse');
        end
    end
pause
end


%% ============= Reliability test 2, partition trials/reps ==========
% change #selected trials, by partitioning 
nPart = [3 5 11];
%%
for i = 1:length(nPart)
    ind_subset = cell(1,nPart(i));
    ind_left = 1:para.nStim;
    for j = 1:nPart(i)
        if j == nPart(i)
            ind_subset{j} = ind_left;
        else
            ind_subset{j} = datasample(ind_left, round(para.nStim/nPart(i)), 'replace',false);
            ind_left = setdiff(ind_left, ind_subset{j});
        end
    end

    
    for k = 1:nPart(i)
        opt.ampLimit    = 0.4.*[0 1];
        opt.trials =    ind_subset{k}; 
%         opt.reps        = Reps_temp(k,:);
%         opt.reps = datasample(1:para.nRep, Reps(i), 'replace',false)
    %     figure, set(gcf, 'color','w')
    %     [X, DataMat_norm] = ViewData(DataMat, para, opt); % X may contain NaNs if there are masked pixels
        [X, DataMat_norm] = getX(DataMat, para, opt);

        % ======== construct a X without NaN ========
        [~, ind_delete]     = find( isnan(X) ); % linear index
        ind_save            = setdiff(1:para.width*para.height, ind_delete);
        X(:, ind_delete)    = [];

        % ======== perform ICA analysis =============
        addpath('D:\=code=\McdermottLab\toolbox_nonparametric-ICA') % MAKE SURE
        % THIS FOLDER IS ON THE TOP OF THE PATH LIST IN MATLAB!!!!
        K = 6;
        RANDOM_INITS = 10;      
        PLOT_FIGURES = 0;
        % ========= reverse the sign for X for intrinsic imaging ========
        tic,[R, W2{i,k}] = nonparametric_ica(X, K, RANDOM_INITS, PLOT_FIGURES);toc
        disp(['i=',num2str(i),', k=',num2str(k)])
    end
end



%% find match with original weights, for partition method
nPart = 3;
P = perms(1:K); % all possible orders
% figure
% set(gcf,'color','w','position', [1323, 380, 1365, 192]);
for m = 1
    for n = 1:nPart
        Dist = zeros(1,size(P,1));
        W_temp = W2{m,n};
        for j = 1:size(P,1)
            ind = P(j,:);
            W_temp2 = W_temp(ind,:);
            Dist(j) = sum(vecnorm(W_temp2' - W_ori'));
        end
        [minDist, indDist] = min(Dist);
        ind = P(indDist,:); % the best matching indecies
        MinDist(m,n) = minDist; % the best matching similarity index (distance)
        % W1 = W(P(indDist,:),:);

        % ======== plot ============
        set(gcf,'color','w','position', [1323, 380, 1365, 192]);
        p = [1 6];
        ha = tight_subplot(p(1),p(2),[0 0],[.1 .01],[.01 .01]);
        for i = 1:K
            cutoff = 0.2;
            comp{i} = zeros(para.height, para.width);
            comp{i}(ind_save) = W_temp(ind(i),:);
            axes(ha(i)); imagesc(comp{i},cutoff.*[-1 1]),axis image, colormap(jet)
            axis off
        end
        pause
    end
end

%% find match with original weights (save reliability measures for each component)
MinDist = cell(1,K);
for i = 1:K
MinDist{i} = zeros(length(Reps), nRep);
end
P = perms(1:K);
% figure
% set(gcf,'color','w','position', [1323, 380, 1365, 192]);

for m = 1:length(Reps)
    for n = 1:nRep
        Dist = zeros(1,size(P,1));
        W_temp = W2{m,n};
        for j = 1:size(P,1)
            ind = P(j,:);
            W_temp2 = W_temp(ind,:);
            Dist(j) = sum(vecnorm(W_temp2' - W_ori'));
        end
        [minDist, indDist] = min(Dist);
        ind = P(indDist,:);
        
        for k = 1:K
            MinDist{k}(m,n) = vecnorm(W_temp(ind(k),:) - W_ori(k,:)) ./ vecnorm(W_ori(k,:));
        end
%         D(m,n) = minDist;
        % W1 = W(P(indDist,:),:);

        % ======== plot ============
        
%         set(gcf,'color','w','position', [1323, 380, 1365, 192]);
%         p = [1 6];
%         ha = tight_subplot(p(1),p(2),[0 0],[.1 .01],[.01 .01]);
%         for k = 1:K
%             cutoff = 0.2;
%             comp{k} = zeros(para.height, para.width);
%             comp{k}(ind_save) = W_temp(ind(k),:);
%             axes(ha(k)); imagesc(comp{k},cutoff.*[-1 1]),axis image, colormap(jet)
%             axis off
%         end
%         pause
    end
end
%% 
figure, set(gcf,'color','w')
for k = 1:K
X = Reps;
Y = mean(MinDist{k},2);
E = std(MinDist{k},[],2)./sqrt(size(MinDist{k},2));
errorbar(X, Y, E,'linewidth',2), hold on
end
xlabel('Percentage of trials selected')
ylabel('Mean Euclidean distance of components')
axis square
legend({'component 1', 'component 2','component 3','component 4','component 5','component 6'})