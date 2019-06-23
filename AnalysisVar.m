function [Var] = AnalysisVar(para,DataMat,opt)
% DataMat = [rep, trial, height, width, frams]
% opt.plotON, opt.ampLimit

% Construct a matrix V to compute variance across REP and STIM
% [rep, stim, pixels(w*h)]
V = zeros(para.nRep,para.nStim,para.width*para.height);
parfor i = 1:para.nRep
    mov =   DataMat(i,:,:,:,:);
    % remove only the 1st dimension (do not use squeeze, in case the stim
    % number is 1 too)
    z =     size(mov);
    mov =   reshape(mov,[z(2:end),1]);
    
    img_base =      squeeze(mean(mov(:,:,:,1:floor(para.fr*para.preStim)),4)); % first second: baseline
    if isempty(opt.tWindow)
        % choose the the response window (from stimulus start to simulus end)
        img_amp =       squeeze(mean(mov(:,:,:,floor(para.fr*para.preStim)+1:floor(para.fr*para.preStim+ + para.fr*para.durStim)),4));
    
    else
        % choose the silence window (prestimulus)
        img_amp =       squeeze(mean(mov(:,:,:,floor(para.fr*opt.tWindow(1))+1:floor(para.fr*opt.tWindow(2))),4));
    end
    img_relamp =    (img_amp - img_base)./img_base;
    
    V(i,:,:) = reshape(img_relamp,para.nStim,para.width*para.height);
end
if para.nRep ~= 1
Var_rep = var(V);
Var_rep = reshape(squeeze(Var_rep),para.nStim,para.height,para.width); % [stim, height, width]
end
if para.nStim ~= 1
Var_stim = var(V,0,2);
Var_stim = reshape(squeeze(Var_stim),para.nRep,para.height,para.width); % [rep, height, width]
end

% Compute variance across frames
Var_frame = zeros(para.nStim,para.height,para.width);
Mean_frame = zeros(para.nStim,para.height,para.width);
Base_frame = zeros(para.nStim,para.height,para.width);
%%
for i = 1:para.nStim
% for i = 1:10
    mov_temp = mean(squeeze(DataMat(:,i,:,:,:)),1); % 1 x height x width x frames
    mov_temp = squeeze(mov_temp); % height x width x frames
    img_base = squeeze(mean(mov_temp(:,:,1:floor(para.fr*para.preStim)),3));
    Base_frame(i,:,:) = img_base; % temporary - keep base frames for each stimulus
    
    img_base = repmat(img_base,1,1,para.nFrame);
    mov_rel = (mov_temp - img_base)./img_base;
%     Var_frame(i,:,:) = var(mov_rel(:,:,floor(para.fr*para.preStim)+1:para.fr*(para.preStim+para.durStim)),0,3);
    Var_frame(i,:,:) = var(mov_rel(:,:,para.fr*opt.tWindow(1)+1:floor(para.fr*opt.tWindow(2))),0,3);
    Mean_frame(i,:,:) = mean(mov_rel(:,:,para.fr*opt.tWindow(1)+1:floor(para.fr*opt.tWindow(2))),3);
    
end

%% see variance map for each stimulus
% for i = 1:10
% figure,imagesc(squeeze(Var_frame(i,:,:))),axis image
% end
%% averaged variance of reps across reps
% Var_rep_mean = mean(Var_rep,1);
% Var_rep_mean = squeeze(Var_rep_mean);
if para.nRep ~= 1
    Var.RepMean =   squeeze(mean(Var_rep,1));
    Var.RepRel =    Var.RepMean./squeeze(mean(img_relamp,1));    % the relative variance (normalized by mean amplitude)

    if opt.plotON
        figure,imagesc(Var.RepMean),axis image, colorbar
        title('Variance across reps, averaged across frames and stimulus')
        figure, imagesc(Var.RepRel,[0 opt.ampLimit]),axis image, colorbar
        title('Variance across reps, averaged across frames and stimulus (normalized)')
    end
end
%% averaged variance of stims across stimuli
% Var_stim_mean = mean(Var_stim,1);
% Var_stim_mean = squeeze(Var_stim_mean);
% figure,imagesc(Var_stim_mean),axis image, colorbar
% title('Variance across stim, averaged across frames and reps')
% % the relative variance (normalized by mean amplitude)
% Var_stim_rel = Var_stim_mean./squeeze(mean(img_relamp,1));
% figure, imagesc(Var_stim_rel,[0 0.1]),axis image, colorbar
% title('Variance across stim, averaged across frames and reps (normalized)')
if para.nStim ~= 1
    Var.StimMean =   squeeze(mean(Var_stim,1));
    Var.StimRel =    Var.StimMean./squeeze(mean(img_relamp,1));    % the relative variance (normalized by mean amplitude)

    if opt.plotON
        figure,imagesc(Var.StimMean),axis image, colorbar
        title('Variance across stimuli, averaged across frames and stimulus')
        figure, imagesc(Var.StimRel,[0 opt.ampLimit]),axis image, colorbar
        title('Variance across stimuli, averaged across frames and stimulus (normalized)')
    end
end
%% averaged variance of frames, across stim and reps
% Var_frame_mean = mean(Var_frame,1);
% Var_frame_mean = squeeze(Var_frame_mean);
% figure,imagesc(Var_frame_mean),axis image, colorbar
% title('Variance across frames, averaged across stimulus and reps')
% % the relative variance (normalized by mean amplitude)
% Var_frame_rel = Var_frame_mean./squeeze(mean(img_relamp,1));
% figure, imagesc(Var_frame_rel,[0 0.06]),axis image, colorbar
% title('Variance across frames, averaged across stimulus and reps (normalized)')

Var.FrameMean =   squeeze(mean(Var_frame,1));
Var.FrameRel =    Var.FrameMean./squeeze(mean(img_relamp,1));    % the relative variance (normalized by mean amplitude)

if opt.plotON
    figure,imagesc(Var.FrameMean),axis image, colorbar
    title('Variance across frames, averaged across frames and stimulus')
    figure, imagesc(Var.FrameRel,[0 opt.ampLimit]),axis image, colorbar
    title('Variance across frames, averaged across frames and stimulus (normalized)')
end
end