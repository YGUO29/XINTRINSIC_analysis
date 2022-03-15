function plotRegressCoef(output, F, para)

Title = {'Frequency (Hz)',...
    'Temp. Mod. (Hz)',...
    'Spec. Mod. (Cyc/Oct)',...
    'Frequency and SpecTemp. Mod.',...
    'Category',...
    'Category (adjusted colorbar)'};
% Title = {'Frequency',...
%     'Frequency and Temp. Mod.',...
%     'Frequency and Spec. Mod.',...
%     'Frequency and SpecTemp. Mod.',...
%     'Category',...
%     'Category (adjusted colorbar)'};
figurex([563         253        1565         883])

for i = [1:3, 6]
    % calculated weighted average of features (the feature with highest
    % coefficients)
    b = output.Coef{i}(:,2:end); % nPix x nFeature
    corr = output.Rsquared(:,i);
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
        case 6 % category
            cb_labels = F.C.category_labels;
            feat = 1:length(cb_labels);
    end
      
    % ==== weighted average as a visualization measurement (negative coefficients problem!!)
%     b = -b;
%     b(b<0) = 0;
%     b = b./repmat(sum(b,2), 1, size(b,2)); 
%     temp = b*[1:length(feat)]'; 
    % ==== max as a visualization measurement ====
    if i == 6
    subplot(2,2,4)
    else
    subplot(2,2,i)
    end
    [~,b_ind] = max(b,[],2); 
    hue = reshape(b_ind,[para.height, para.width]);
    
    % set saturation to one for all pixels 
%     corr = ones(1,para.height*para.width);
    
    % find pixels that have all zero coefficients
    ind_zero = find(sum(b,2) == 0);
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