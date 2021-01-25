function myButtonDownFcn(src, ~)
    global H D
    H.cp = get(gca, 'CurrentPoint');
    H.cp = round(H.cp);
    H.cp
    
    D.img_withMark = D.img;
    D.img_withMark(H.cp(1,2)-1:H.cp(1,2)+1, H.cp(1,1)-1:H.cp(1,1)+1) = NaN;
    set(H.img, 'CData', D.img_withMark)
    
    
D.variable.trace = squeeze(D.DataMat_norm(:,H.cp(1,2), H.cp(1,1), :))';
plot(H.panel.trace, D.variable.trace);
title(H.panel.trace, 'Temporal trace')

D.variable.tuning_2d = squeeze(D.TuningMap(:,:,sub2ind([D.para.height, D.para.width],  H.cp(1,2), H.cp(1,1))));
set(H.tuning_2d, 'CData', D.variable.tuning_2d);

D.variable.tuning_F1 = mean(squeeze(D.TuningMap(:,:,sub2ind([D.para.height, D.para.width],  H.cp(1,2), H.cp(1,1)))), 2);
cla(H.panel.tuning_marginal1)
semilogx(H.panel.tuning_marginal1, D.F1, D.variable.tuning_F1);
title(H.panel.tuning_marginal1, 'Frequency tuning (marginal)')

D.variable.tuning_F2 = mean(squeeze(D.TuningMap(:,:,sub2ind([D.para.height, D.para.width],  H.cp(1,2), H.cp(1,1)))), 1);
cla(H.panel.tuning_marginal2)
plot(H.panel.tuning_marginal2, fliplr(D.variable.tuning_F2));
title(H.panel.tuning_marginal2, 'Level tuning (marginal)')


    
    %% plot temporal trace
%     trace = squeeze(D.DataMat_norm(:, H.cp(1,2), H.cp(1,1), :)); % flip x and y
%     plot(H.plot_trace, trace');
%     set(H.plot_trace, 'YLim', [-0.05, 0.25])
%     title('Frequency-Level tuning map')
% 
%     
%     % plot tuning map
%     set(H.plot_tuning, 'CData', squeeze(D.TuningMap(:,:,sub2ind([D.para.height, D.para.width], H.cp(1,2), H.cp(1,1)))) )
%     axis image
%     
%     % feature 1/2 tuning
%     tuning_F1 = mean(squeeze(D.TuningMap(:,:,sub2ind([D.para.height, D.para.width],  H.cp(1,2), H.cp(1,1)))), 2);
%     semilogx(H.plot_tuning_marginal1, D.F1, tuning_F1);
%     title('Frequency tuning (marginal)')
% 
%     tuning_F2 = mean(squeeze(D.TuningMap(:,:,sub2ind([D.para.height, D.para.width],  H.cp(1,2), H.cp(1,1)))), 1);
%     plot(H.plot_tuning_marginal2, D.F2, tuning_F2);
%     title('Level tuning (marginal)')

    
    
    
end