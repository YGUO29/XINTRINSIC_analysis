function plotRegressCoef2(output, F, para)

%% plot coefficients from case #5 (combined freq+mod features)
Title = {'Frequency (Hz)',...
    'Temp. Mod. (Hz)',...
    'Spec. Mod. (Cyc/Oct)'};
figurex([563         253        1565         883])
idx = cell(1,3);
idx{1} = 1:F.nFreq;
idx{2} = 1+F.nFreq : F.nFreq+F.nTemp;
idx{3} = 1+F.nFreq+F.nTemp : F.nFreq+F.nTemp+F.nSpec;
for i = 1:3
    % calculated weighted average of features (the feature with highest
    % coefficients)
    b = output.Coef{5}(:,2:end); % nPix x nFeature
    corr = output.Rsquared(:,5);
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
    end
      
    % ==== weighted average as a visualization measurement (negative coefficients problem!!)
%     b = -b;
%     b(b<0) = 0;
%     b = b./repmat(sum(b,2), 1, size(b,2)); 
%     temp = b*[1:length(feat)]'; 
    % ==== max as a visualization measurement ====
    subplot(1,3,i)
    [~,b_ind] = max(b(:,idx{i}),[],2); 
    hue = reshape(b_ind,[para.height, para.width]);
    
    % set saturation to one for all pixels 
%     corr = ones(1,para.height*para.width);
    
    % find pixels that have all zero coefficients
    ind_zero = find(sum(b(:,idx{i}),2) == 0);
    corr(ind_zero) = 0;
    
    sat = reshape(corr,[para.height, para.width]);
    para_plot.continuous_cb = 0;
    para_plot.hue_max       = length(feat);
    para_plot.mirror        = 0;
    para_plot.cb_labels     = cb_labels;
    para_plot.title         = Title{i};
    plotMap_with2Dcolorbar(sat, hue, para_plot)

    if isfield(para, 'ct')
        plotContour(para.ct);
    end
    
    if para.mirror
        set(gca,'XDir','reverse')
    end
end