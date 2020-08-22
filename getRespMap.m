function X = getRespMap(DataMat_norm, para, options)
if ~isfield(para, 'frames')
    para.frames = 1:para.nFrame;
end
    
if strcmp(options.method, 'mean')
    % ========= use mean amplitude as a measurement ====
    X = reshape(mean(X1,1), para.height, para.width);
    figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w');
    imagesc(X), axis image, colormap(jet)
    title(options.title)
elseif strcmp(options.method, 'var')
    % ==== use variance as a measurement ====
    temp = permute(DataMat_norm,[2 3 4 1]);
    temp = temp(:,:,para.frames,:);
    temp = reshape(temp,[para.height, para.width, size(DataMat_norm,1)*length(para.frames)]);
    X = var(temp,[],3);
    figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w');
    if isfield(options,'X_exp')
        Max = max(log10(options.X_exp(:)));Min = min(log10(options.X_exp(:)));
    else
        Max = max(log10(X(:)));Min = min(log10(X(:)));
    end
    imagesc(log10(X), [Min Max]), colormap(jet), axis image; 
    h = colorbar; title(h,'          10^')
    title(options.title)
else
    disp('set 3rd variable as mean or var')
end
end