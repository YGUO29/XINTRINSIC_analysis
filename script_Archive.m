%% perform ICA analysis (Norman-Haignere 2015)
addpath('D:\=code=\McdermottLab\toolbox_nonparametric-ICA') % MAKE SURE
% THIS FOLDER IS ON THE TOP OF THE PATH LIST IN MATLAB!!!!
K = 6;
RANDOM_INITS = 10;      
PLOT_FIGURES = 0;
% ========= reverse the sign for X for intrinsic imaging ========
tic,[R, W] = nonparametric_ica(-X, K, RANDOM_INITS, PLOT_FIGURES);toc
%% NMF 
K = 20;
X_pos = -X - min(-X(:));
opt = statset('Maxiter',200);
[R, W] = nnmf(X_pos, K, 'options', opt);
%% run multiple K values
addpath('D:\=code=\McdermottLab\toolbox_nonparametric-ICA') % MAKE SURE
% THIS FOLDER IS ON THE TOP OF THE PATH LIST IN MATLAB!!!!
nK = 20;
X_hat = cell(1, nK);
for K = 1:nK
    RANDOM_INITS = 10;      
    PLOT_FIGURES = 0;
    % ========= reverse the sign for X for intrinsic imaging ========
    tic,[R, W] = nonparametric_ica(X, K, RANDOM_INITS, PLOT_FIGURES);toc
    X_hat{K} = R * W;
end

% calculate MSE for the entire X (all stims, all pixels)
error = zeros(1, nK); 
for K = 1:nK
    
    error(K) = norm(X - X_hat{K}, 'fro');
    
end
figurex;
plot(mean(error, 1),'marker','*')
xlabel('# components')
ylabel('reconstruction error')

% calculate MSE for each row of X (for each stim, all pixels)
% error = zeros(para.nStim, nK); 
% for K = 1:nK
%     for i = 1:para.nStim
%         error(i, K) = immse(X(i,:), X_hat{K}(i,:));
%     end
% end
% figurex;
% errorbar(mean(error, 1), std(error, [], 1),'linewidth',2)
% xlabel('# components')
% ylabel('MSE')

%% tight plot
ind = 1:K;
% % ind = [3 2 4 6 5 1]; % for 80Z
% ind = [1 2 5 4 3 6]; % for 132D 1/19 session
% ind = [1 2 4 6 5 3]; % for 132D, session 2
p = [1, K];
figurex;
set(gcf,'color','w','position', [1323, 380, 1365, 192]);
ha = tight_subplot(p(1),p(2),[.01 .01],[.1 .01],[.01 .01]);
for i = 1:K
    cutoff = 0.08;
    cutoff = mean(W(i,:)) + 7*std(W(i,:)); % variable cutoff values for each components
    comp{i} = zeros(para.height, para.width);
    comp{i}(ind_save) = W(ind(i),:);
% %     % ============= plot components only =============
%     mask_temp = double(~mask_outline_reg);
%     mask_temp(mask_temp == 0) = -inf;
    axes(ha(i)); imagesc(comp{i},cutoff.*[-1 1]),axis image, colormap(jet)
    colorbar,
    axis off
%     plotContour(ct)
end