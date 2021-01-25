function [R, W, comp, recon_error, X_hat] = runDecomp(X, Ks, opt, para)
% opt.fluo = 0; 
% opt.method = 'mICA'; % or 'NMF'
% opt.nRows = 1;
% opt.plotON = 1;
% Ks = 6;
% [Rs, Ws, comps, recon_error, X_hat] = runDecomp(X, Ks, opt, para);

% initialization 
nK = length(Ks);
R = cell(1, nK);
W = R;
comp = R;
X_hat = cell(1, nK);
recon_error = zeros(1, nK);

% reverse X if it's intrinsic (although this doesn't seem to change the
% results for mICA)
if ~opt.fluo 
    X = -X;
    if isfield(opt, 'X_test'); opt.X_test = -opt.X_test; end
else
end

% figurex;
% ha = tight_subplot(max(Ks), max(Ks), [.01 .01],[.1 .01],[.01 .01]);
if strcmp(opt.method, 'PCA')
    k = 1; % store Rs and Ws in only one cell
    [U, S, V] = svd(X, 'econ');
    R{k} = U(:, 1:Ks)';
    W{k} = V(:, 1:Ks)';
    recon_error = diag(S)./sum(diag(S));
    if opt.plotON
        % for plot components.
        figurex([1323, 380, 1365, 192]); 
        ha = tight_subplot(opt.nRows, ceil(Ks./opt.nRows), [.01 .01],[.1 .01],[.01 .01]);

        % plot components
        for i = 1:Ks
%             cutoff = 0.02;
            cutoff = mean(W{k}(i,:)) + 7*std(W{k}(i,:)); % variable cutoff values for each components
            comp_temp = zeros(para.height, para.width);
            comp_temp(para.ind_save) = W{k}(i,:);
            % ============= plot components only =============
            %     mask_temp = double(~mask_outline_reg);
            %     mask_temp(mask_temp == 0) = -inf;
            axes(ha(i)); 
%             axes(ha(i+(k-1)*max(Ks)));
            imagesc(comp_temp,cutoff.*[-1 1]),axis image, colormap(jet)
            comp{k}(:, :, i) = comp_temp;
            drawnow;
            colorbar;
            axis off
            %     plotContour(ct)
        end
    end
else
    
%% run a loop over nK    
for k = 1:nK
    K = Ks(k); % k: the index; K: the #components   
    switch opt.method 
        case 'mICA' % modified ICA (Sam, 2015)
            addpath('D:\SynologyDrive\=code=\McdermottLab\toolbox_nonparametric-ICA') % MAKE SURE
            % THIS FOLDER IS ON THE TOP OF THE PATH LIST IN MATLAB!!!!
            RANDOM_INITS = 10;      
            PLOT_FIGURES = 0;
            % ========= reverse the sign for X for intrinsic imaging ========
            tic,[R{k}, W{k}] = nonparametric_ica(X, K, RANDOM_INITS, PLOT_FIGURES);toc
            
            % adjust order of W by absolute values
            [~, ind] = sort( mean(abs(W{k}), 2), 'descend');
            W{k} = W{k}(ind, :);
            R{k} = R{k}(:, ind);
%             X_hat{k} = R{k} * W{k};
            if isfield(opt, 'X_test')
                X_test = opt.X_test;
            else
                X_test = X;
            end
            X_hat{k} = (R{k} / (R{k}' * R{k}))* R{k}' * X_test;
            recon_error(k) = norm(X_test - X_hat{k}, 'fro');
        case 'NMF'
            X_pos = X - min(X(:));
            options = statset('Maxiter',200);
            [R{k}, W{k}] = nnmf(X_pos, K, 'options', options);
%             X_hat{k} = R{k} * W{k};
            if isfield(opt, 'X_test')
                X_test = opt.X_test - min(opt.X_test(:));
            else
                X_test = X_pos;
            end
            X_hat{k} = (R{k} / (R{k}' * R{k}))* R{k}' * X_test;
            recon_error(k) = norm(X_test - X_hat{k}, 'fro');
    end
    
    
    if opt.plotON
        % for plot components.
        figurex; 
        set(gcf,'color','w','position', [1323, 380, 1365, 192]);
        ha = tight_subplot(opt.nRows, ceil(K./opt.nRows), [.01 .01],[.1 .01],[.01 .01]);

        % arrange components' order
        ind = 1:K;
        % % ind = [3 2 4 6 5 1]; % for 80Z
        % ind = [1 2 5 4 3 6]; % for 132D 1/19 session
        % ind = [1 2 4 6 5 3]; % for 132D, session 2

        % plot components
        for i = 1:K
            cutoff = 0.02;
            cutoff = mean(W{k}(ind(i),:)) + 7*std(W{k}(ind(i),:)); % variable cutoff values for each components
            comp_temp = zeros(para.height, para.width);
            comp_temp(para.ind_save) = W{k}(ind(i),:);
            % ============= plot components only =============
            %     mask_temp = double(~mask_outline_reg);
            %     mask_temp(mask_temp == 0) = -inf;
            axes(ha(i)); 
%             axes(ha(i+(k-1)*max(Ks)));
            imagesc(comp_temp,cutoff.*[-1 1]),axis image, colormap(jet)
            comp{k}(:, :, i) = comp_temp;
            drawnow;
            colorbar;
            axis off
            %     plotContour(ct)
        end
    end
end
end

% plot reconstruction error
figurex;
plot(Ks, recon_error,'marker','*')
xlabel('# components')
ylabel('reconstruction error')

end