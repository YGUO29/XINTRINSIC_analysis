function matched_comp = runMatchSubjects(Decomp1, Decomp2, opt, para)

Ws_all = [];
Rs_all = [];
for i = 1:length(opt.Ks)
    Ws_all = cat(1, Ws_all, Decomp2.Ws{i});
    Rs_all = cat(2, Rs_all, Decomp2.Rs{i});
end

% find the minimal distance match
% [comp_dist, comp_ind] = pdist2(Rs_all', Decomp1.Rs{1}', 'euclidean', 'Smallest',1);

% find the maximum correlation match
for i = 1: size(Decomp1.Rs{1},2)
temp = corrcoef([Decomp1.Rs{1}(:,i), Rs_all]);
[comp_dist(i), comp_ind(i)] = max(temp(1,2:end));
end

% convert index in the concatenated matrix to [i, K] (ith component when
% component number = K)
% assuming chosen K = 5, 6, 7 
matched_comp = zeros(3, size(Decomp1.Rs{1},2));
for i = 1:size(Decomp1.Rs{1},2)
    ind_Ks = min( find(comp_ind(i) - cumsum(opt.Ks) - 1<0) );
    matched_comp(1,i) = ind_Ks;
    matched_comp(2,i) = opt.Ks(ind_Ks);
    if ind_Ks == 1
        matched_comp(3,i) = comp_ind(i);
    else
        matched_comp(3,i) = comp_ind(i) - sum(opt.Ks(1:ind_Ks-1));
    end
end

figurex([1440        1171         560         167]); 
for i = 1:size(Decomp1.Rs{1},2)
    subplot(1,2,1)
    imagesc(Decomp1.comp{1}(:,:,i)), axis image, axis off, colormap(jet)
    set(gca,'XDir','reverse')
    subplot(1,2,2)
%     plotTonotopyWithComponents(1, Decomp_80Z, T, para)
    imagesc(reshape((Decomp2.Ws{matched_comp(1,i)}(matched_comp(3,i),:)),...
        para.height, para.width)), axis image, axis off, colormap(jet)
    pause
end