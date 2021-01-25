% Methods to determine the number of components, by calculating variance
% explained by these components
opt = struct;
% opt.tWindow = [para.preStim, para.preStim + para.durStim];
if ~isfield(opt, 'tWindow')
    opt.tWindow = [para.preStim, para.preStim + para.durStim];
end

set_V = nchoosek(1:para.nRep,2);
% set_V = [5 6];
% ============= changed permutation number here ========
nPerm = size(set_V,1);
% =====================================================

result = zeros(para.nRep, 2, nPerm); % second dimension: measure 1 and measure 2 
for K = 1:para.nRep % number of components
    for iPerm = 1:nPerm
        
        % use the rest of the reps to calculate R and W
        opt.reps = setdiff(1:para.nRep, set_V(iPerm,:));
%         opt.reps = datasample(1:para.nRep, para.nRep-2, 'replace', false);
%         set_V = setdiff(1:para.nRep, opt.reps);
%         ind_V = datasample(1:5, 2, 'replace', false);

        X = getX(DataMat, para, opt);
        RANDOM_INITS = 10;      
        PLOT_FIGURES = 0;
        [~, ind_delete]     = find( isnan(X) ); % linear index
        ind_save            = setdiff(1:para.width*para.height, ind_delete);
        X(:, ind_delete)    = [];
        % === ====== reverse the sign for X for intrinsic imaging ========
        tic,[R, W] = nonparametric_ica( X, K, RANDOM_INITS, PLOT_FIGURES);toc

        V = zeros(2, para.nStim, para.width*para.height);
%         batch{1} = datasample(1:para.nRep, 5, 'replace', false);
%         batch{2} = setdiff(1:para.nRep, batch{1});
        ind_V = set_V(iPerm,:);
        
        for i = 1:2
            ind = ind_V(i);
            mov = squeeze(DataMat(ind,:,:,:,:));
            img_base = squeeze(mean(mov(:,:,:,1:floor(para.fr*para.preStim)),4)); % first second: baseline
            img_amp = squeeze(mean(mov(:,:,:,floor(para.fr*opt.tWindow(1))+1 : floor(para.fr*opt.tWindow(2))),4));
            % figure,imshow(img_base,[])
            % figure,imshow(img_amp,[])
            img_relamp = (img_amp - img_base)./img_base;
            V(i,:,:) = reshape(img_relamp,para.nStim,para.width*para.height);
        end



        V1 = squeeze(V(1,:,:)); 
        V2 = squeeze(V(2,:,:)); 

        rho = zeros(para.width*para.height,1);
        rho_norm = rho;
        for i = 1:para.height*para.width
            v1 = V1(:,i); v2 = V2(:,i);
            v_proj1 = R*pinv(R)*v1;
            v_proj2 = R*pinv(R)*v2;
            rho1 = corr(v_proj1, v2);
            rho2 = corr(v_proj2, v1);
            rho(i) = tanh(0.5*(atanh(rho1)+atanh(rho2)));

            rho_norm(i) = rho(i)/sqrt(corr(v1,v2)*corr(v_proj1,v_proj2));
            rho_norm(i) = abs(rho_norm(i)).^2;
        %     rho(i) = corr(v_proj1,v_proj2).^2;
        end
        result(K,1,iPerm) = median(rho);
        result(K,2,iPerm) = median(rho_norm);
        iPerm, K, 

    end
   
end

%%
figure, set(gcf,'color','w')
Y = mean(result,3);
E = std(result,[],3);
errorbar(Y(:,1),E(:,1)./sqrt(nPerm),'linewidth',2,'color','k');ylabel('Prediction Accuracy','fontsize',14); xlabel('# Components','fontsize',14)
set(gca,'fontsize',20), axis square
figure, set(gcf,'color','w')
errorbar(Y(:,2),E(:,2)./sqrt(nPerm),'linewidth',2,'color','k');ylabel('% of Explained Variance','fontsize',14); xlabel('# Components','fontsize',14)
set(gca,'fontsize',20), axis square
%%
% Methods to determine the number of components, by calculating variance
% explained by these components
opt = struct;
% opt.tWindow = [para.preStim, para.preStim + para.durStim];
if ~isfield(opt, 'tWindow')
    opt.tWindow = [para.preStim, para.preStim + para.durStim];
end

set_V = nchoosek(1:para.nRep,2);
% set_V = [5 6];
nPerm = size(set_V,1)/2;
result = zeros(para.nRep, 2, nPerm); % second dimension: measure 1 and measure 2 
for K = 1:para.nRep % number of components
    for iPerm = 1:nPerm
        
        % use the rest of the reps to calculate R and W
        opt.reps = setdiff(1:para.nRep, set_V(iPerm,:));
%         opt.reps = datasample(1:para.nRep, para.nRep-2, 'replace', false);
%         set_V = setdiff(1:para.nRep, opt.reps);
%         ind_V = datasample(1:5, 2, 'replace', false);

        X = getX(DataMat, para, opt);
        RANDOM_INITS = 10;      
        PLOT_FIGURES = 0;
        [~, ind_delete]     = find( isnan(X) ); % linear index
        ind_save            = setdiff(1:para.width*para.height, ind_delete);
        X(:, ind_delete)    = [];
        % === ====== reverse the sign for X for intrinsic imaging ========
        tic,[R, W] = nonparametric_ica( X, K, RANDOM_INITS, PLOT_FIGURES);toc

        V = zeros(2, para.nStim, para.width*para.height);
%         batch{1} = datasample(1:para.nRep, 5, 'replace', false);
%         batch{2} = setdiff(1:para.nRep, batch{1});
        ind_V = set_V(iPerm,:);
        
        for i = 1:2
            ind = ind_V(i);
            mov = squeeze(DataMat(ind,:,:,:,:));
            img_base = squeeze(mean(mov(:,:,:,1:floor(para.fr*para.preStim)),4)); % first second: baseline
            img_amp = squeeze(mean(mov(:,:,:,floor(para.fr*opt.tWindow(1))+1 : floor(para.fr*opt.tWindow(2))),4));
            % figure,imshow(img_base,[])
            % figure,imshow(img_amp,[])
            img_relamp = (img_amp - img_base)./img_base;
            V(i,:,:) = reshape(img_relamp,para.nStim,para.width*para.height);
        end


        V1 = squeeze(V(1,:,:)); 
        V2 = squeeze(V(2,:,:)); 

        rho = zeros(para.width*para.height,1);
        rho_norm = rho;
        
        V1_proj =  R*pinv(R)*V1;
        V2_proj =  R*pinv(R)*V2;
        
        tic 
        for i = 1:para.height*para.width
            v1 = V1(:,i); v2 = V2(:,i);
            v1_proj = V1_proj(:,i); v2_proj = V2_proj(:,i);
             rho1(i) = corr(v1_proj, v2);
             rho2(i) = corr(v2_proj, v1);
        end
        toc
       
        rho_mean = tanh(0.5*(atanh(rho1)+atanh(rho2)));

        rho_norm = rho_mean./sqrt(corr(V1, V2).*corr(V1_proj, V2_proj));

        result(K,1,iPerm) = median(rho);
        result(K,2,iPerm) = median(rho_norm);
        iPerm, K, 

    end
   
end