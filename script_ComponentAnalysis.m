% Component analysis including reliability analysis
%% View data, get the original X
opt = struct;
% opt.trials = 28:36; 
[X, DataMat_norm] = getX(DataMat, para, opt);

opt.ampLimit    = 0.4.*[0 1];
%     figure, set(gcf, 'color','w')
%     [X, DataMat_norm] = ViewData(DataMat, para, opt); % X may contain NaNs if there are masked pixels
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
tic,[R, W] = nonparametric_ica(X, K, RANDOM_INITS, PLOT_FIGURES);toc
ind = [3 2 4 6 5 1]; % for 80Z
% ind = [1 2 5 4 3 6]; % for 132D 1/19 session
% ind = [1 2 4 6 5 3]; % for 132D, session 2
W_ori = W(ind,:);

%% ============= Reliability test 1, select Perc of the trials/reps ==========
% change #selected trials, by randomly select certain percent
Perc = [fliplr(0.7:0.05:0.95), 0.5, 0.3, 0.1, 0.035];
nRep = 10;
% W1 = cell(length(Perc), nRep);

% change number of reps
Reps = 1:4;
W2 = cell(length(Reps), nRep);
%%
% Reps_temp = zeros(nRep, 4);
% for i = 1:nRep
%   Reps_temp(i,:) =  datasample(1:para.nRep, 4, 'replace',false);
% end

for i = length(Reps)
% for i = 1:length(Perc)
    for k = 1:nRep
        opt.ampLimit    = 0.4.*[0 1];
%         opt.trials =    datasample(1:para.nStim,round(Perc(i)*para.nStim),'replace',false); 
        opt.reps        = Reps_temp(k,:);
%         opt.reps = datasample(1:para.nRep, Reps(i), 'replace',false)
    %     opt.tWindow     = [para.preStim, para.preStim + 18]; % start and end of integration window for calculating response amplitude
        opt.tWindow     = []; % start and end of integration window for calculating response amplitude
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
        tic,[R, W2{i,k}] = nonparametric_ica(-X, K, RANDOM_INITS, PLOT_FIGURES);toc
        disp(['Perc=',num2str(i),', Rep=',num2str(k),' Selected ',num2str(round(Perc(i)*para.nStim)),' trials'])

    end
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
%% tight plot
comp = cell(1,K);
% ind = 1:6;
ind = [3 2 4 6 5 1]; % for 80Z
% ind = [1 2 5 4 3 6]; % for 132D 1/19 session
% ind = [1 2 4 6 5 3]; % for 132D, session 2figure

set(gcf,'color','w','position', [1323, 380, 1365, 192]);
% ha = tight_subplot(p(1),p(2),[.01 .01],[.1 .01],[.01 .01]);
p = [1 6];
ha = tight_subplot(p(1),p(2),[0 0],[.1 .01],[.01 .01]);
for i = 1:K
    cutoff = 0.2;
%     cutoff = mean(W(i,:)) + 7*std(W(i,:)); % variable cutoff values for each components
    comp{i} = zeros(para.height, para.width);
    comp{i}(ind_save) = W(ind(i),:);
% %     % ============= plot components only =============
    axes(ha(i)); imagesc(comp{i},cutoff.*[-1 1]),axis image, colormap(jet)
%     colorbar,
    axis off
end

%% find match with original weights
D = zeros(length(Reps), nRep);
P = perms(1:K); % all possible orders
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
        D(m,n) = minDist;
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
        D(m,n) = minDist; % the best matching similarity index (distance)
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
%% 
X = Reps;
Y = mean(D,2);
E = std(D,[],2);
figure, set(gcf,'color','w')
errorbar(X, Y, E,'linewidth',2)
xlabel('Percentage of trials selected')
ylabel('Mean Euclidean distance of components')
axis square

%% find match with original weights (save reliability measures for each component)
D = cell(1,K);
for i = 1:K
D{i} = zeros(length(Reps), nRep);
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
            D{k}(m,n) = vecnorm(W_temp(ind(k),:) - W_ori(k,:)) ./ vecnorm(W_ori(k,:));
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
Y = mean(D{k},2);
E = std(D{k},[],2)./sqrt(size(D{k},2));
errorbar(X, Y, E,'linewidth',2), hold on
end
xlabel('Percentage of trials selected')
ylabel('Mean Euclidean distance of components')
axis square
legend({'component 1', 'component 2','component 3','component 4','component 5','component 6'})