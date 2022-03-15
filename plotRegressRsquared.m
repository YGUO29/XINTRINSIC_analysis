function plotRegressRsquared(output, F, para)

%% R-squared map
% Corr1 = output.Rsquared(:,1:6);
Corr1 = output.Rsquared(:,1:6);
Corr2 = output.Rsquared(:,6);
Title = {'Frequency',...
    'Temp. Mod.',...
    'Spec. Mod.',...
    'Frequency + SpecTemp. Mod.',...
    'Frequency + Spec + Temp. Mod.',...
    'Category',...
    'Category (adjusted colorbar)'};
% Title = {'Frequency',...
%     'Frequency and Temp. Mod.',...
%     'Frequency and Spec. Mod.',...
%     'Frequency and SpecTemp. Mod.',...
%     'Category',...
%     'Category (adjusted colorbar)'};
figurex([563         253        1565         553])
Max = max(Corr1(:));
% Max = 0.3;
Min = prctile(Corr1(:),0);
for i = 1:6
    subplot(2,3,i)
    b_ind = reshape(Corr1(:,i),[para.height, para.width]);
%     Min = prctile(temp(:),75);
    Min = 0;
    
    imagesc(b_ind, [Min, Max]), axis image, colorbar
    title(Title{i})
    axis off
    colormap(jet)
    if isfield(para, 'ct')
        plotContour(para.ct, 'RdBu');
    end
    
    if para.mirror
        set(gca,'XDir','reverse')
    end
end

figurex([1031         711        1756         627])
Max = max(output.Rsquared(:));
Min = prctile(Corr2(:),0);
for i = 1:2
    subplot(1,2,i)
    if i== 2
        b_ind = reshape(output.Rsquared(:,6),[para.height, para.width]);
        Max = max(b_ind(:));
%         Min = prctile(Corr(:),0);
        Min = 0;
    else
        b_ind = reshape(output.Rsquared(:,6),[para.height, para.width]);
    end
%     Max = max(temp(:));
    
    imagesc(b_ind, [Min, Max]), axis image, colorbar
    title(Title{i+5})
    axis off
    colormap(jet)
    
    if isfield(para, 'ct')
        plotContour(para.ct, 'RdBu');
    end
    
    if para.mirror
        set(gca,'XDir','reverse')
    end
end

% plot R^2 map with tonotopy
figurex([1031         711        1756         627])
Max = max(output.Rsquared(:));
Min = prctile(Corr2(:),0);
for i = 1:2
    subplot(1,2,i)
    if i== 2
        b_ind = reshape(output.Rsquared(:,6),[para.height, para.width]);
        Max = max(b_ind(:));
%         Min = prctile(Corr(:),0);
        Min = 0;
    else
        b_ind = reshape(output.Rsquared(:,6),[para.height, para.width]);
    end
%     Max = max(temp(:));

    imagesc(b_ind, [Min, Max]), axis image, colorbar
    title(Title{i+5})
    axis off
    colormap(jet)
    
    if isfield(para, 'ct')
        plotContour(para.ct, 'RdBu');
    end
    
    if para.mirror
        set(gca,'XDir','reverse')
    end
end