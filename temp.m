Max = max(abs(DataMat_mean(:)));
figure, imagesc(DataMat_mean, [-Max, Max]), axis image, colormap(jet), colorbar
h = images.roi.Circle(gca,'Center',[floor(para.width/2) floor(para.height/2)],'Radius',floor(para.height/2)); 
title('Press Enter after the position is adjusted')
pause
mask = createMask(h);
%%
figure,
DataMat_test1 = squeeze(DataMat_norm(1,:,:,3:end)).*repmat(mask,1,1,348);
plot(mean(reshape(DataMat_test1,75*120, 348), 1)), hold on

DataMat_test2 = squeeze(DataMat_norm(2,:,:,1:end-2)).*repmat(mask,1,1,348);
plot(mean(reshape(DataMat_test2,75*120, 348), 1))

%% get sound features and category features
nGroup = max(F.table.cat_number);
coch_env_group = cell(1,nGroup);
coch_env_mean = zeros(size(F.coch_env, 1), nGroup);
for i = 1:size(F.table,1)
    ind = F.table.cat_number(i);
    coch_env_group{ind} = [coch_env_group{ind}, F.coch_env(:,i)];
end
figurex;
for i = 1:nGroup
    coch_env_mean(:, i) = mean(coch_env_group{i}, 2);
    semilogx(F.cf_log,  coch_env_mean(:, i), 'color', F.C.colors(i,:)), hold on
end
legend(F.C.category_labels)
%% get average coch_env
coch_env_mean = zeros(size(F.coch_env, 1), 1);
for i = 1:15
    coch_env_mean = coch_env_mean + F.coch_env(:,ind_descend(i));
end
coch_env_mean = coch_env_mean ./ 15;
figurex;
semilogx(F.cf_log,  coch_env_mean./max(coch_env_mean)), hold on
semilogx(F.cf_log,  mean(F_voc.coch_env, 2)./max(mean(F_voc.coch_env, 2)))


    
%% calculate weight of cochlear envelope
coch_weight = zeros(1, size(F.table,1));
for i = 1:size(F.table,1)
    coch_weight(i) = (F.cf_log*F.coch_env(:,i)) ./ sum(F.coch_env(:,i));
end

[~, ind_descend] = sort(coch_weight, 'descend');
figurex; 
hist(coch_weight)
    