% analysis for Ripple sounds (trial based)

%% Calculate MTF matrix (nSM x nTM x nPix)
S.sm            = [0, 1/8, 1/4, 1/2, 1, 2, 4, 8];
S.tm            = [1/2, 1, 2, 4, 8, 16, 32, 64];
S.tm            = [fliplr(-S.tm), 0, S.tm]; % all TM rates
nSM = length(S.sm);
nTM = length(S.tm);
nPix = size(X,2);
MTF     = zeros(nSM,nTM,nPix);


TM      = S.tm; nTM = length(TM);
SM      = S.sm; nSM = length(SM);
IndMat  = cell(nSM, nTM); 
for i = 1:nSM
    for j = 1:nTM
        IndMat{i,j} = [SM(i), TM(j)];
    end
end
%     IndPix  = find(mask_final); % linear index of selected ROI
%     IndPix  = find(idx == kk);
%     IndPix = 1:size(X,2);
for kk = 1:nPix
    order   = [nTM:-2:1,2:2:nTM-1]; % rearrange the order to align with sounds (in the experiments sound TM order is 0, 1/2, -1/2, ...)
    data_temp   = X(:,kk); % select a pixel. rows of X: number of SM/TM parameters
    data_temp   = - reshape(data_temp, nTM, nSM)'; % each column --> same SM --> transpose
    MTF(:,:,kk)  = data_temp;
    MTF(:,:,kk)  = MTF(:,order,kk);

end
MTF_mean = mean(MTF,3);

%% show response patterns in a MTF grid
figurex;
cmax = prctile(MTF(:),99);
cmin = prctile(MTF(:),1);
for iSM = 1:nSM
    for iTM = 1:nTM
        subplot(nSM, nTM, nTM*(iSM-1)+iTM)
        imagesc(reshape(MTF(iSM,iTM,:), para.height, para.width), [cmin, cmax])
        colormap(jet)
        axis off
        axis image
    end
end
        
%% get best TM and SM map
bTM = zeros(1, nPix);
bSM = bTM;
for i = 1:nPix
    MTF_temp = MTF(:,:,i);
    [row, col] = find(MTF_temp == max(MTF_temp(:)));
%     bSM(i) = SM(row);
%     bTM(i) = TM(col);
    bSM(i) = row;
    bTM(i) = col;

end
bSM = reshape(bSM, para.height, para.width);
bTM = reshape(bTM, para.height, para.width);

figurex;
imagesc(bSM), colormap(jet(8)), axis image, colorbar
title('Best spectral modulation')

figurex;
imagesc(bTM), colormap(jet(17)), axis image, colorbar
title('Best temporal modulation')
%% plot all MTF for all pixels
figurex([1440         351        1066         987]); 
imagesc(flipud(MTF_mean)); colormap(jet), colorbar, axis image
set(gca,'XTickLabel',arrayfun(@num2str,TM(2:2:16),'UniformOutput',false))
set(gca,'yTickLabel',arrayfun(@num2str,fliplr(SM(2:2:8)),'UniformOutput',false))
xlabel('Temporal Modulation Rate (Hz)')
ylabel('Spectral Modulation Rate (Cyc/Oct)')
pause

mm = prctile(MTF(:),100);
h = imagesc(flipud(MTF(:,:,1)),[0 mm]); colormap(jet), colorbar, axis image
set(gca,'XTickLabel',arrayfun(@num2str,TM(2:2:16),'UniformOutput',false))
set(gca,'yTickLabel',arrayfun(@num2str,fliplr(SM(2:2:8)),'UniformOutput',false))
xlabel('Temporal Modulation Rate (Hz)')
ylabel('Spectral Modulation Rate (Cyc/Oct)')
for i = 1:nPix
    % plot averaged MTF over pixels
    pause(0.01);
    set(h,'CData',flipud(MTF(:,:,i)));
end

%% plot component's response profile as MTF 
% rearrange Rs to get MTF of each component
K = 10;
R = Rs{1};
figurex([1440         351        1066         987]); 
for i = 1:K
    subplot(2,4,i)
    mm = prctile(R(:),95);
    MTF_temp = R(:,i); % select a pixel. rows of X: number of SM/TM parameters
    MTF_temp = reshape(MTF_temp, nTM, nSM)';
    order   = [nTM:-2:1,2:2:nTM-1];
    MTF_temp  = MTF_temp(:,order);
    
    h = imagesc(flipud(MTF_temp),[0 mm]); colormap(parula), colorbar, axis image
    xticks( 2:7:16 )
    xticklabels(arrayfun(@num2str,TM(2:7:16),'UniformOutput',false))
    yticks(2:2:8)
    yticklabels(arrayfun(@num2str,fliplr(SM(1:2:8)),'UniformOutput',false))
%     xlabel('Temporal Modulation Rate (Hz)')
%     ylabel('Spectral Modulation Rate (Cyc/Oct)')
end

%% manually circle a region, plot MTF
I = comp{1}(:,:,3);
figure,imshow(I,[]);
h = impoly(gca)
mask_final = createMask(h);
IndPix  = find(mask_final); % linear index of selected ROI
MTF_temp = mean(MTF(:,:,IndPix), 3);

figurex([1440         351        1066         987]); 

mm = prctile(MTF_temp(:),100);
h = imagesc(flipud(MTF_temp),[0 mm]); colormap(parula), colorbar, axis image
xticks( 2:7:16 )
xticklabels(arrayfun(@num2str,TM(2:7:16),'UniformOutput',false))
yticks(1:2:8)
yticklabels(arrayfun(@num2str,fliplr(SM(1:2:8)),'UniformOutput',false))
xlabel('Temporal Modulation Rate (Hz)')
ylabel('Spectral Modulation Rate (Cyc/Oct)')


%% classify MTFs
k = 10; % k-means cluster of MTF
idx = kmeans(X',k);
I = rescale(I,0,1);
f1 = figure;
for i = 1:k
    ind = find(idx == i);

    II = repmat(I,1,1,3); % RGB plot for pixel positions
    Itemp = ones(para.height, para.width);
    Itemp(ind) = 0;
    II(:,:,1) = II(:,:,1) .* Itemp;
    figure(f1),subplot(1,k,i),imagesc(II),axis image
    
end