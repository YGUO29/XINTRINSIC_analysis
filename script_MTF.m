figure,imshow(I,[]);
h = impoly(gca)
mask_final = createMask(h);
%% script for calculating MTF
S.fm            = [0, 1/8, 1/4, 1/2, 1, 2, 4, 8];
S.tm            = [1/2, 1, 2, 4, 8, 16, 32, 64];
S.tm            = [fliplr(-S.tm), 0, S.tm]; % all frequencies
% for kk = 1:5
    TM      = S.tm; nTM = length(TM);
    SM      = S.fm; nSM = length(SM);
    IndMat  = cell(nSM, nTM); 
    for i = 1:nSM
        for j = 1:nTM
            IndMat{i,j} = [SM(i), TM(j)];
        end
    end


    IndPix  = find(mask_final); % linear index of selected ROI
%     IndPix  = find(idx == kk);
%     IndPix = 1:size(X,2);
    nPix    = length(IndPix);
    MTF     = zeros(nSM,nTM,nPix);
    order   = [nTM:-2:1,2:2:nTM-1];
    for i = 1:nPix
        data_temp   = X(:,IndPix(i)); % rows of X: number of SM/TM parameters
        MTF(:,:,i)  = - reshape(data_temp, nSM, nTM);
        MTF(:,:,i)  = MTF(:,order,i);
    end
    MTF_mean = mean(MTF,3);

    % plot averaged MTF over pixels
    mm = prctile(MTF_mean(:),100);
    figure,
    h = imagesc(flipud(MTF_mean),[0 mm]); colormap(hot), colorbar, axis image
    set(gca,'XTickLabel',arrayfun(@num2str,TM(2:2:16),'UniformOutput',false))
    set(gca,'yTickLabel',arrayfun(@num2str,fliplr(SM(1:2:7)),'UniformOutput',false))
    xlabel('Temporal Modulation Rate (Hz)')
    ylabel('Frequency Modulation Rate (Cyc/Oct)')
% end
%%
pause 
for i = 1:nPix
    set(h,'CData',flipud(MTF(:,:,i)));
    pause(0.1)
end

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