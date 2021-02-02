% view data interactively
global H D

% D.DataMat_norm = DataMat_norm; % contains temporal info for all pixels
% D.RMap = RMap; % #Feature1 x #Feature2 x #pixels
% D.para = para;
% D.F1 = F1; D.F2 = F2;
% D.F1_log = log2(D.F1./110)+2; % log scale, number means A-x (starting from A2)
%%
H.fig = figurex([17         274        1895         825]);

% the image to click on
D.img = reshape(mean(X, 1), D.para.height, D.para.width);

subplot(2,3,1)
Max = max(abs(D.img(:)));
H.img = imagesc(D.img, [-Max, Max]); axis image
set(H.img, 'ButtonDownFcn', @myButtonDownFcn)
title('Mean response amplitude')

% ==== temporal traci ====
H.panel.trace = subplot(2,3,[2,3]);
D.variable.trace = squeeze(D.DataMat_norm(:, D.para.height/2, D.para.width/2, :))';
H.trace = plot(H.panel.trace, D.variable.trace);
set(H.panel.trace, 'YLim', [-0.05, 0.25])
title('Temporal trace')

% ==== 2D tuning map ====
H.panel.tuning_2d = subplot(2,3,4);
D.variable.tuning_2d = squeeze(D.RMap(:,:,sub2ind([D.para.height, D.para.width],  D.para.height/2, D.para.width/2)));
H.tuning_2d = imagesc(flipud(D.variable.tuning_2d), 'XData', D.F2, 'YData', D.F1);
% H.tuning_2d = imagesc(D.variable.tuning_2d, 'XData', D.F2, 'YData', D.F1_log);
colorbar
title('2D tuning map')
xlabel('Sound level')
ylabel('Frequency (units = A-x)')

% ==== tuning for F1 ====
H.panel.tuning_marginal1 = subplot(2,3,5);
D.variable.tuning_F1 = mean(squeeze(D.RMap(:,:,sub2ind([D.para.height, D.para.width],  D.para.height/2, D.para.width/2))), 1);
if isfield(D,'F1_log')
    H.tuning_marginal1 = plot(D.F1_log, D.variable.tuning_F1);
else
    H.tuning_marginal1 = plot(D.variable.tuning_F1);
end
xticks(D.ticks{1}); xticklabels(D.ticklabels{1})
title('Marginal tuning for Feature #1')

% ==== tuning for F2 ====
H.panel.tuning_marginal2 = subplot(2,3,6);
D.variable.tuning_F2 = mean(squeeze(D.RMap(:,:,sub2ind([D.para.height, D.para.width],  D.para.height/2, D.para.width/2))), 2);
H.tuning_marginal2 = plot(D.F2, D.variable.tuning_F2);
title('Marginal tuning for Feature #2')



