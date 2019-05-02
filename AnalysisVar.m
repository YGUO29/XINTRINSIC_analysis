function [Var_rep_mean, Var_rep_rel, Var_stim_mean, Var_stim_rel] = AnalysisVar(para,DataMat)

V = zeros(para.nRep,para.nStim,para.width*para.height);
for i = 1:para.nRep
    mov = squeeze(DataMat(i,:,:,:,:));

    img_base = squeeze(mean(mov(:,:,:,1:floor(para.fr*para.preStim)),4)); % first second: baseline
    img_amp = squeeze(mean(mov(:,:,:,floor(para.fr*para.preStim)+1:floor(para.fr*para.preStim+ + para.fr*para.durStim)),4));
    % figure,imshow(img_base,[])
    % figure,imshow(img_amp,[])
    img_relamp = (img_amp - img_base)./img_base;
    V(i,:,:) = reshape(img_relamp,para.nStim,para.width*para.height);
end

Var_rep = var(V);
Var_rep = reshape(squeeze(Var_rep),para.nStim,para.height,para.width);
Var_stim = var(V,0,2);
Var_stim = reshape(squeeze(Var_stim),para.nRep,para.height,para.width);

% compute variance across frames
Var_frame = zeros(para.nStim,para.height,para.width);

for i = 1:para.nStim
    mov_temp = mean(squeeze(DataMat(:,i,:,:,:)),1); % height x width x frames
    mov_temp = squeeze(mov_temp);
    img_base = squeeze(mean(mov_temp(:,:,1:floor(para.fr*para.preStim)),3));
    img_base = repmat(img_base,1,1,para.nFrame);
    mov_rel= (mov_temp - img_base)./img_base;
    Var_frame(i,:,:) = var(mov_rel(:,:,floor(para.fr*para.preStim)+1:para.fr*(para.preStim+para.durStim)),0,3);
end
%% see variance map for each stimulus
% for i = 1:10
% figure,imagesc(squeeze(Var_frame(i,:,:))),axis image
% end
%% averaged variance of reps across stimulus
Var_rep_mean = mean(Var_rep,1);
Var_rep_mean = squeeze(Var_rep_mean);
figure,imagesc(Var_rep_mean),axis image, colorbar
title('Variance across reps, averaged across frames and stimulus')
% the relative variance (normalized by mean amplitude)
Var_rep_rel = Var_rep_mean./squeeze(mean(img_relamp,1));
figure, imagesc(Var_rep_rel,[0 0.1]),axis image, colorbar
title('Variance across reps, averaged across frames and stimulus (normalized)')
%% averaged variance of stims across reps
Var_stim_mean = mean(Var_stim,1);
Var_stim_mean = squeeze(Var_stim_mean);
figure,imagesc(Var_stim_mean),axis image, colorbar
title('Variance across stim, averaged across frames and reps')
% the relative variance (normalized by mean amplitude)
Var_stim_rel = Var_stim_mean./squeeze(mean(img_relamp,1));
figure, imagesc(Var_stim_rel,[0 0.1]),axis image, colorbar
title('Variance across stim, averaged across frames and reps (normalized)')

%% averaged variance of frames, across stim and reps
Var_frame_mean = mean(Var_frame,1);
Var_frame_mean = squeeze(Var_frame_mean);
figure,imagesc(Var_frame_mean),axis image, colorbar
title('Variance across frames, averaged across stimulus and reps')
% the relative variance (normalized by mean amplitude)
Var_frame_rel = Var_frame_mean./squeeze(mean(img_relamp,1));
figure, imagesc(Var_frame_rel,[0 0.06]),axis image, colorbar
title('Variance across frames, averaged across stimulus and reps (normalized)')
end