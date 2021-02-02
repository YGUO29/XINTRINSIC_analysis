
%% Parameters for frequency-level tuning
global D

octaves = 0:1/12:(8+1/4); % one semitone apart
freqs = 110.*2.^octaves;
% Feature 1 changes first
D.F1 = freqs;
D.F2 = 0:10:80;
D.F2 = fliplr(D.F2);
D.fluo = 1;

D.Titles{'Best Frequencies', 'Best Sound Level'};
D.ticks{1} = 0:8;
D.ticklabels{1} = {'A2', 'A3', 'A4', 'A5', 'A6', 'A7', 'A8', 'A9', 'A10'};
% plot_opt.ticks{2} = 0:8;
% plot_opt.ticklabels{2} = {'A2', 'A3', 'A4', 'A5', 'A6', 'A7', 'A8', 'A9', 'A10'};
%% Parameters for MTF
% F1 = TM, which changes first
D.F1            = [1/2, 1, 2, 4, 8, 16, 32, 64, 128]; 
% D.F1          = [1/2, 1, 2, 4, 8, 16, 32, 64]; 
D.F1          = [fliplr(-D.F1), 0, D.F1]; % all TM rates
D.F1_order    = [length(D.F1):-2:1,2:2:length(D.F1)-1]; 
D.F2          = [0, 1/8, 1/4, 1/2, 1, 2, 4, 8];
D.F2_order    = 1:length(D.F2);
D.fluo        = 1;

D.Titles = {'Best temporal modulation', 'Best spectral modulation'};
D.ticks{1} = 1:2:length(S.tm);
D.ticklabels{1} = arrayfun(@num2str,D.F1(1:2:end),'UniformOutput',false);
D.ticks{2} = 1:length(S.sm);
D.ticklabels{2} = arrayfun(@num2str,D.F2,'UniformOutput',false);


%%
D.nF1 = length(D.F1);
D.nF2 = length(D.F2);
nPix = size(X,2);
D.RMap = zeros(D.nF2, D.nF1, nPix); % 2nd feature is on the 1st dimension 

%     IndPix  = find(mask_final); % linear index of selected ROI
%     IndPix  = find(idx == kk);
%     IndPix = 1:size(X,2);
for kk = 1:nPix
    if D.fluo % positive response
        data_temp   = X(:,kk); % select a pixel. rows of X: number of SM/TM parameters
    else % negative response
        data_temp   = - X(:,kk); % select a pixel. rows of X: number of SM/TM parameters
    end
    data_temp   = reshape(data_temp, D.nF1, D.nF2)'; % each column --> same SM --> transpose
    
    data_temp   = data_temp(:, D.F1_order); % permute to have the same order as F1 labels
    data_temp   = data_temp(D.F2_order, :); % permute to have the same order as F2 labels
    
    D.RMap(:,:,kk)  = data_temp;
end
RMap_mean = mean(D.RMap,3);
if D.fluo
    X_maxproj = reshape(max(X, [], 1), para.height, para.width);
else
     X_maxproj = reshape(max(-X, [], 1), para.height, para.width);
end
%% get best F1 and F2
bF1 = zeros(1, nPix);
bF2 = bF1;
for i = 1:nPix
    RMap_temp = D.RMap(:,:,i);
    % winner takes all
    [row, col] = find(RMap_temp == max(RMap_temp(:)));    
    % marginal + winner takes all
%     mF1 = mean(RMap_temp,1); % marginal tuning for F1
%     [~, col] = max(mF1);
%     mF2 = mean(RMap_temp,2); % marginal tuning for F2
%     [~, row] = max(mF2);

    if length(row) > 1 || length(col) > 1
        bF1(i) = NaN;
        bF2(i) = NaN;
    else
        bF1(i) = col;
        bF2(i) = row;
    end

    % marginal + weighted average
%     mF1 = mean(RMap_temp,1); % marginal tuning for F1
%     mF1(mF1<0) = 0;
%     weights = mF1./sum(mF1);
% %     bF1(i) = F1*weights;
%     bF1(i) = (1:D.nF1)*weights';
%     
%     mF2 = mean(RMap_temp,2); % marginal tuning for F2
%     mF2(mF2<0) = 0;
%     weights = mF2./sum(mF2);
% %     bF2(i) = F2*weights';
%     bF2(i) = (1:D.nF2)*weights;
end
bF1 = reshape(bF1, para.height, para.width);
bF2 = reshape(bF2, para.height, para.width);
%% plot bF1 and bF2 maps (without amplitude weighting)
figurex;
imagesc(bF1), axis image, colorbar
% imagesc(log2(bF1./110)), axis image, colorbar
% colormap(jet(D.nF1)), 
% colormap(parula(D.nF1))
CT = cbrewer('div','RdBu',D.nF1); colormap(CT)
caxis([0 D.nF1])
colorbar('Ticks',D.ticks{1},...
         'TickLabels',D.ticklabels{1});
title(D.Titles{1})

% ==== unify positive and negative TMs ====
figurex;
bF1(bF1<9) = 18 - bF1((bF1<9));
imagesc(bF1), axis image, colorbar
% imagesc(log2(bF1./110)), axis image, colorbar
% colormap(jet(D.nF1)), 
colormap(parula((D.nF1+1)/2))
caxis([floor(D.nF1/2) D.nF1])
colorbar('Ticks',D.ticks{1},...
         'TickLabels',D.ticklabels{1});
title(D.Titles{1})

figurex;
imagesc(bF2), axis image, colorbar
% colormap(jet)
colormap(parula(D.nF2)), 
caxis([0,D.nF2]);
colorbar('Ticks',D.ticks{2},...
         'TickLabels',D.ticklabels{2});
title(D.Titles{2})
%% plot bF1 and bF2 maps (with amplitude weighting)
bF1(isnan(bF1)) = 1;
bF2(isnan(bF2)) = 1;

figurex;
CT = cbrewer('div','RdBu',D.nF1); 
img = CT(bF1,:); img = reshape(img, para.height, para.width, 3);
img = img.*repmat(rescale(X_maxproj, 0.2, 1), 1, 1, 3);

imagesc(img), axis image, colormap(CT), colorbar
caxis([0 D.nF1])
colorbar('Ticks',D.ticks{1},...
         'TickLabels',D.ticklabels{1});
title(D.Titles{1})

% ==== unify positive and negative TMs ====
figurex;
bF1(bF1<D.nF1/2) = D.nF1 + 1 - bF1((bF1<D.nF1/2));
bF1 = bF1 - floor(D.nF1/2);
CT = colormap(parula((D.nF1+1)/2));
img = CT(bF1,:); img = reshape(img, para.height, para.width, 3);
img = img.*repmat(rescale(X_maxproj, 0.2, 1), 1, 1, 3);

imagesc(img), axis image, colormap(CT), colorbar
caxis([floor(D.nF1/2) D.nF1])
colorbar('Ticks',D.ticks{1},...
         'TickLabels',D.ticklabels{1});
title(D.Titles{1})

figurex;
CT = (parula(D.nF2));
img = CT(bF2,:); img = reshape(img, para.height, para.width, 3);
img = img.*repmat(rescale(X_maxproj, 0.2, 1), 1, 1, 3);
imagesc(img), axis image, colormap(CT), colorbar
caxis([0,D.nF2]);
colorbar('Ticks',D.ticks{2},...
         'TickLabels',D.ticklabels{2});
title(D.Titles{2})


%% fit tuning curves with Gaussian
D.fit_para{1} = zeros(nPix, 2); % mean and std
D.fit_para{2} = zeros(nPix, 2); % mean and std

% x1 = [ceil(D.nF1/2):D.nF1]'; 
x1 = [1:ceil(D.nF1/2)]'; 
x2 = [1:D.nF2]';
for i = 1:nPix
    D.variable.tuning_F1 = mean(D.RMap(:,:,i),1);
    y1 = zeros(ceil(D.nF1/2),1);
    y1(1) = D.variable.tuning_F1(ceil(D.nF1/2)); 
    y1(2:end) = mean([D.variable.tuning_F1(1:floor(D.nF1/2)); D.variable.tuning_F1(ceil(D.nF1/2)+1:end)])';
    [f1, gof] = fit(x1, y1, 'gauss1');
    D.fit_para{1}(i,1) = f1.a1;
    D.fit_para{1}(i,2) = f1.b1;
    D.fit_para{1}(i,3) = f1.c1;
    D.fit_para{1}(i,4) = gof.rsquare;
%     figure, plot(f1, x1+floor(D.nF1/2), y1);

    D.variable.tuning_F2 = mean(D.RMap(:,:,i),2);
    y2 = D.variable.tuning_F2;
    [f2, gof] = fit(x2, y2, 'gauss1');
    D.fit_para{2}(i,1) = f2.a1;
    D.fit_para{2}(i,2) = f2.b1;
    D.fit_para{2}(i,3) = f2.c1;
    D.fit_para{2}(i,4) = gof.rsquare;

    if ~mod(i,100)
        i
    end
%     figure, plot(f2, x2, y2);
%     pause
end
%%
figurex;
img = D.fit_para{2}(:,4);
% lm = 0; img(img<lm) = lm; 
% um = 8; img(img>um) = um;
imagesc(reshape(img, para.height, para.width)), axis image
colorbar

%% classify MTFs
I = reshape(mean(X, 1), para.height, para.width);
k = 10; % k-means cluster of MTF
ik = kmeans(X',k); % cluter numbers 1~k
I = rescale(I,0,1);
f1 = figurex;
% f2 = figure;
for i = 1:k
    ind = find(ik == i);

    II = repmat(I,1,1,3); % RGB plot for pixel positions
    Itemp = ones(para.height, para.width);
    Itemp(ind) = 0;
    II(:,:,1) = II(:,:,1) .* Itemp;
    figure(f1),subplot(2,k,i),imagesc(II),axis image
    
    figure(f1), subplot(2,k,i+k)
    MTF_temp = mean(MTF(:,:,ind), 3);
    mm = prctile(MTF_temp(:),100);
    h = imagesc(flipud(MTF_temp),[0 mm]); colormap(parula), colorbar, axis image
    xticks( 2:2:18 )
    xticklabels(arrayfun(@num2str,TM(2:2:18),'UniformOutput',false))
    yticks(1:2:8)
    yticklabels(arrayfun(@num2str,fliplr(SM(1:2:8)),'UniformOutput',false))
    xlabel('Temporal Modulation Rate (Hz)')
    ylabel('Spectral Modulation Rate (Cyc/Oct)')
end