Max = max(abs(DataMat_mean(:)));
figure, imagesc(DataMat_mean, [-Max, Max]), axis image, colormap(jet), colorbar
h = images.roi.Circle(gca,'Center',[floor(para.width/2) floor(para.height/2)],'Radius',floor(para.height/2)); 
title('Press Enter after the position is adjusted')
pause
mask = createMask(h);


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

    
%% calculate weight of cochlear envelope
coch_weight = zeros(1, size(F.table,1));
for i = 1:size(F.table,1)
    coch_weight(i) = (F.cf_log*F.coch_env(:,i)) ./ sum(F.coch_env(:,i));
end

[~, ind_descend] = sort(coch_weight, 'descend');
figurex; 
hist(coch_weight)
    
%% get average coch_env
coch_env_top15 = zeros(size(F.coch_env, 1), 1);
coch_env_voc = mean(F.coch_env(:,find(F.table.cat_number == 12)), 2);
ct = 0;
i = 1;
while ct <= 15
    if F.table.cat_number(ind_descend(i)) ~= 12
        coch_env_top15 = coch_env_top15 + F.coch_env(:,ind_descend(i));
        ct = ct + 1;
    else
    end
    i = i + 1;
end
coch_env_top15 = coch_env_top15 ./ 15;
figurex;
semilogx(F.cf_log,  coch_env_voc), hold on
semilogx(F.cf_log,  coch_env_top15), hold on

legend({'15 voc', '15 Nat w/ highest center of mass'})
