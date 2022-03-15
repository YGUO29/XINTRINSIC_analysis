function plotComp(k, Decomp, opt, para)
K = size(Decomp.Ws{k},1); % number of components
figurex;
set(gcf,'color','w','position', [1323, 380, 1365, 192]);
ha = tight_subplot(opt.nRows, ceil(K./opt.nRows), [.01 .01],[.1 .01],[.01 .01]);

% plot components
ind = 1:K;
for i = 1:K
%             cutoff = 0.15;
    cutoff = mean(Decomp.Ws{k}(ind(i),:)) + 7*std(Decomp.Ws{k}(ind(i),:)); % variable cutoff values for each components
    comp_temp = zeros(para.height, para.width);
    comp_temp(para.ind_save) = Decomp.Ws{k}(ind(i),:);
    % ============= plot components only =============
    %     mask_temp = double(~mask_outline_reg);
    %     mask_temp(mask_temp == 0) = -inf;

    axes(ha(i)); 
%             axes(ha(i+(k-1)*max(opt.Ks)));
    imagesc(comp_temp, cutoff.*[-1 1]),axis image, colormap(jet)
    comp{k}(:, :, i) = comp_temp;
    drawnow;
    colorbar;
    axis off
    if para.mirror
        set(gca,'XDir','reverse');
    else
    end

%     if isfield(para,'ct')
%         plotContour(para.ct, 'gray')
%     end
end