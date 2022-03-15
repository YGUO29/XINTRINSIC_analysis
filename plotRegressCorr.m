function plotRegressCorr(output, F, T, para, animal)
%% correlation between freq and best modulations
% amplitude mask
para = getTonotopyParameters(para, animal);
octaves = 0:1/para.freq_per_oct:para.freq_oct_range; % one semitone apart
freqs = para.freq_base.*2.^octaves; 

mask1_threshold = 0.02;
mask1_ind = find(T.sat_map > mask1_threshold);

% x: bf from tonotopy; y: coefficients
y = zeros(para.height*para.width,3);
x = reshape(T.hue_map, 1, para.height*para.width);
for i = 1:3
    [~,y(:,i)] = max(output.Coef{i}(:,2:end),[],2); 
end

% scatter plot 1, tonotopy and other coefficients and
figurex([432         889        2825         449]);
for i = 1:3
mask2_ind = find(output.Rsquared(:,i) > 0.3);
mask_ind = intersect(mask1_ind, mask2_ind);

subplot(1,3,i)
scatter(x(mask_ind), y(mask_ind,i))
corr = corrcoef(x(mask_ind), y(mask_ind,i));
xticks(0:para.freq_per_oct:length(freqs)-1)
xlim([0 length(freqs)-1])
set(gca,'xticklabels',arrayfun(@num2str, freqs(1:para.freq_per_oct:length(freqs)),'UniformOutput',false))
xlabel('Best frequency (Hz)')
title(['correlation = ', num2str(corr(1,2), '%.2f')])
switch i
    case 1
        ylim([1 F.nFreq])
        ylabel('Coefficients - frequencies')
        yticks(1:F.nFreq)
        set(gca,'yticklabels',arrayfun(@num2str,F.FreqBounds(2:end),'UniformOutput',false))
    case 2
        ylim([1 F.nSpec])
        ylabel('Coefficients - spec. modulation')
        yticks(1:F.nSpec)
        set(gca,'yticklabels',arrayfun(@num2str,F.spec_mod_rates,'UniformOutput',false))
    case 3
        ylim([1 F.nTemp])
        ylabel('Coefficients - temp. modulation')
        yticks(1:F.nTemp)
        set(gca,'yticklabels',arrayfun(@num2str,F.temp_mod_rates(2:end),'UniformOutput',false))
    otherwise
end
end


% modulation trade-off
figurex([1440         821         599         517]); 
histogram2(y(mask_ind,2), y(mask_ind,3));
set(gca,'ytick', 1:F.nSpec, 'YLim', [0.5 F.nSpec+0.5],...
    'yticklabels',arrayfun(@num2str,F.spec_mod_rates,'UniformOutput',false))
% xlabel('Spectral modulation rate (Cyc/Oct)')
set(gca, 'xtick', 1:F.nTemp,'XLim', [0.5 F.nTemp+0.5],...
    'xticklabels',arrayfun(@num2str,F.temp_mod_rates(2:end),'UniformOutput',false))
% ylabel('Temporal modulation rate (Hz)')

% scatter plot 2, tonotopy and other coefficients and
idx = cell(1,3);
idx{1} = 1:F.nFreq;
idx{2} = 1+F.nFreq : F.nFreq+F.nTemp;
idx{3} = 1+F.nFreq+F.nTemp : F.nFreq+F.nTemp+F.nSpec;
b = output.Coef{5}(:,2:end);
figurex([432         889        2825         449]);
for i = 1:3
    [~,y] = max(b(:,idx{i}),[],2); 
    
    mask2_ind = find(output.Rsquared(:,5) > 0.3);
    mask_ind = intersect(mask1_ind, mask2_ind);

    subplot(1,3,i)
    scatter(x(mask_ind), y(mask_ind))
    corr = corrcoef(x(mask_ind), y(mask_ind));
    xticks(0:para.freq_per_oct:length(freqs)-1)
    xlim([0 length(freqs)-1])
    set(gca,'xticklabels',arrayfun(@num2str, freqs(1:para.freq_per_oct:length(freqs)),'UniformOutput',false))
    xlabel('Best frequency (Hz)')
    title(['correlation = ', num2str(corr(1,2), '%.2f')])
    switch i
        case 1
            ylim([1 F.nFreq])
            ylabel('Coefficients - frequencies')
            yticks(1:F.nFreq)
            set(gca,'yticklabels',arrayfun(@num2str,F.FreqBounds(2:end),'UniformOutput',false))
        case 2
            ylim([1 F.nSpec])
            ylabel('Coefficients - spec. modulation')
            yticks(1:F.nSpec)
            set(gca,'yticklabels',arrayfun(@num2str,F.spec_mod_rates,'UniformOutput',false))
        case 3
            ylim([1 F.nTemp])
            ylabel('Coefficients - temp. modulation')
            yticks(1:F.nTemp)
            set(gca,'yticklabels',arrayfun(@num2str,F.temp_mod_rates(2:end),'UniformOutput',false))
        otherwise
    end
end



% R-squared mask
figurex([644         857        1356         481]);
for i = 1:3
ind = reshape(output.Rsquared(:,1:3),[para.height, para.width, 3]);
[~, ind_max] = max(ind,[],3);
mask2_ind = find(ind_max == i);

mask_ind = intersect(mask1_ind, mask2_ind);
test = zeros(para.height, para.width);
test(mask_ind) = 1;

subplot(3,2,i*2), imshow(test)
if para.mirror
   set(gca,'XDir','reverse') 
end

subplot(3,2,i*2-1), hist(T.hue_map(mask_ind),20)
xlim([0 length(freqs)-1])
xticks(0:para.freq_per_oct:length(freqs)-1)
set(gca,'xticklabels',arrayfun(@num2str, freqs(1:para.freq_per_oct:length(freqs)),'UniformOutput',false))
if i == 3
    xlabel('Best frequency (Hz)')
end

end

end
