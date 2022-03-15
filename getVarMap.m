function Var = getVarMap(para,DataMat,opt)
% compute variance across repetitions, stimulus and frames

% output: a structure containing absolute variance and relative variance in
% these conditions

% usage:
% opt = struct;
% opt.ampLimit    = [0 98]; % the percentile of the data distribution as display limits
% opt.tWindow     = [4 5]; % or don't set this up, for default tWindow
% Var = getVarMap(para,DataMat,opt)
% DataMat is a matrix of [rep, stim, height, width, frams]

if ~isfield(opt,'tWindow')
    opt.tWindow = [para.preStim para.durStim + para.preStim];
end
if ~isfield(opt,'plotON')
    opt.plotON = 1;
end

%% Var_rep and Var_stim: based on response pattern
% average relative response amplitude within tWindow 
% generate a matrix of [rep, stim, pixels(w*h)]
V = zeros(para.nRep, para.nStim, para.width*para.height);
for i = 1:para.nRep
    mov     = DataMat(i,:,:,:,:);
    % remove only the 1st dimension 
    % (do not use squeeze, in case the stim number is 1 too)
    z       = size(mov);
    mov     = reshape(mov,[z(2:end),1]);
    
    img_base    = squeeze(mean(mov(:,:,:,1:floor(para.fr*para.preStim)),4)); % first second: baseline
    % choose the the response window (from stimulus start to simulus end)
    img_amp     = squeeze(mean(mov(:,:,:,floor(para.fr*opt.tWindow(1))+1:floor(para.fr*opt.tWindow(2))),4));
    img_relamp  = (img_amp - img_base)./img_base;
    V(i,:,:)    = reshape(img_relamp,para.nStim,para.width*para.height);
end

% variance across repetitions
if para.nRep ~= 1
    Var_rep = var(V);
    Var_rep = reshape(squeeze(Var_rep),para.nStim,para.height,para.width); % [stim, height, width]
end

% variance across stimulus
if para.nStim ~= 1
    Var_stim = var(V,0,2);
    Var_stim = reshape(squeeze(Var_stim),para.nRep,para.height,para.width); % [rep, height, width]
% end
%% method 2: Var_rep and Var_stim: based on response traces
% (did not perform well, it seems to be very sensitive to noise)
% For each sound and each pixel, take the temporal traces from all repetitions 
% (within stimulus-presentation window) and calculate a similarity among these traces

% mov_rel_sep = zeros(size(DataMat));
% for i = 1:para.nRep
%     for j = 1:para.nStim
%     img_base =  squeeze(mean(DataMat(i, j, :, :, 1:floor(para.fr*para.preStim)), 5, 'omitnan')); % first second: baseline
%     % choose the the response window (from stimulus start to simulus end)
% %     img_amp     = squeeze(mean(mov(:,:,:,floor(para.fr*opt.tWindow(1))+1:floor(para.fr*opt.tWindow(2))),4));
%     img_base = repmat(img_base, 1, 1, para.nFrame);
%     mov_rel_sep(i,j,:,:,:) = (squeeze((DataMat(i,j,:,:,:))) ...
%         - img_base)./img_base;     
%     end
% end
% mov_rel_sep = mov_rel_sep(:,:,:,:,floor(para.fr*opt.tWindow(1))+1:floor(para.fr*opt.tWindow(2)));
% img_relamp = squeeze(mean(mov_rel_sep,5));
% img_relamp = squeeze(mean(img_relamp,1));
% 
% V = reshape(mov_rel_sep, para.nRep, para.nStim, para.width*para.height, ...
%     size(mov_rel_sep,5));
% 
% for iPix = 1:para.height*para.width
%     for iStim = 1:para.nStim
%         Var_rep(iStim, iPix) = mean(pdist(squeeze(V(:,iStim,iPix,:)))); % nRep x nFrames
%     end
% end
% Var_rep = reshape(squeeze(Var_rep),para.nStim,para.height,para.width);
% 
% for iPix = 1:para.height*para.width
%     for iRep = 1:para.nRep
%         Var_stim(iRep, iPix) = mean(pdist(squeeze(V(iRep,:,iPix,:)))); % nRep x nFrames
%     end
% end
% Var_stim = reshape(squeeze(Var_stim),para.nRep,para.height,para.width);

%% method 1: Compute variance across frames
Var_frame = zeros(para.nStim,para.height,para.width);
Mean_frame = zeros(para.nStim,para.height,para.width);
% Base_frame = zeros(para.nStim,para.height,para.width); % temporary

for i = 1:para.nStim
% for i = 1:10
    % pick a stimulus, average across rep
    mov_temp = mean(squeeze(DataMat(:,i,:,:,:)),1); % 1 x height x width x frames
    mov_temp = squeeze(mov_temp); % height x width x frames
    img_base = squeeze(mean(mov_temp(:,:,1:floor(para.fr*para.preStim)),3));
%     Base_frame(i,:,:) = img_base; % temporary - keep base frames for each stimulus
    
    img_base = repmat(img_base,1,1,para.nFrame);
    mov_rel = (mov_temp - img_base)./img_base;
    Var_frame(i,:,:) = var(mov_rel(:,:,para.fr*opt.tWindow(1)+1:floor(para.fr*opt.tWindow(2))),0,3);
    Mean_frame(i,:,:) = mean(mov_rel(:,:,para.fr*opt.tWindow(1)+1:floor(para.fr*opt.tWindow(2))),3);
end

Mean_frame = squeeze(mean(Mean_frame,1)); 
Mean_frame(Mean_frame<0) = NaN;
%% method 2: Compute variance across frames
% Var_frame = zeros(para.nRep, para.nStim, para.height,para.width);
% Mean_frame = zeros(para.nRep, para.nStim,para.height,para.width);
% % Base_frame = zeros(para.nStim,para.height,para.width); % temporary
% 
% for iStim = 1:para.nStim
%     for iRep = 1:para.nRep
% % for i = 1:10
%     % pick a stimulus,and a  rep
%     mov_temp = squeeze(DataMat(iRep, iStim,:,:,:)); % 1 x height x width x frames
%     img_base = squeeze(mean(mov_temp(:,:,1:floor(para.fr*para.preStim)),3));
% %     Base_frame(i,:,:) = img_base; % temporary - keep base frames for each stimulus
%     
%     img_base = repmat(img_base,1,1,para.nFrame);
%     mov_rel = (mov_temp - img_base)./img_base;
%     Var_frame(iRep, iStim, :,:) = var(mov_rel(:,:,para.fr*opt.tWindow(1)+1:floor(para.fr*opt.tWindow(2))),0,3);
%     Mean_frame(iRep, iStim, :,:) = mean(mov_rel(:,:,para.fr*opt.tWindow(1)+1:floor(para.fr*opt.tWindow(2))),3);
%     end
% end
% Var_frame = squeeze(mean(Var_frame,1));
% Mean_frame = squeeze(mean(Mean_frame,[1,2])); Mean_frame(Mean_frame<0) = NaN;
%% see variance map for each stimulus
% for i = 1:10
% figure,imagesc(squeeze(Var_frame(i,:,:))),axis image
% end
fig1 = figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w','position',[943 600 2304 420]);
fig2 = figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w','position',[943 600 2304 420]);

%% averaged variance of reps across reps
if para.nRep ~= 1
    % mean variance across stimulus
    Var.RepMean = squeeze(mean(Var_rep,1));
    % the mean response amplitude across stimulus
%     temp = abs(squeeze(mean(img_relamp,1))); temp(temp<0) = NaN;

    % relative variance, convert to log scale for visualization
    Var.RepRel = log10(Var.RepMean./Mean_frame);    % the relative variance (normalized by mean amplitude)
    Var.RepMean = log10(Var.RepMean);

    if opt.plotON
        figure(fig1), subplot(1,3,1)
        if opt.flip
            imagesc(fliplr(Var.RepMean)),axis image, 
        else
            imagesc(Var.RepMean),axis image, 
        end
        h = colorbar; title(h,'          10^')
        caxis([prctile(Var.RepMean(:),opt.ampLimit(1)) prctile(Var.RepMean(:),opt.ampLimit(2))])
        axis off
        title('Variance across repetitions \newline (averaged across stimulus)','fontsize',16)
        
        figure(fig2), subplot(1,3,1)
        if opt.flip
            imagesc(fliplr(Var.RepRel)),axis image, 
        else
            imagesc(Var.RepRel),axis image, 
        end
        h = colorbar; title(h,'          10^')
        caxis([prctile(Var.RepRel(:),opt.ampLimit(1)) prctile(Var.RepRel(:),opt.ampLimit(2))])
        axis off
        title('Normalized variance across repetitionss \newline (averaged across stimulus)', 'fontsize',16)
    end
end
%% averaged variance of stims across stimuli
if para.nStim ~= 1
    Var.StimMean =   squeeze(mean(Var_stim,1));
%     temp = abs(squeeze(mean(img_relamp,1))); temp(temp<0) = NaN;
    Var.StimRel =    log10(Var.StimMean./Mean_frame);    % the relative variance (normalized by mean amplitude)
    Var.StimMean = log10(Var.StimMean);
    
    if opt.plotON
        figure(fig1), subplot(1,3,2), 
        if opt.flip
            imagesc(fliplr(Var.StimMean)),axis image, 
        else
            imagesc(Var.StimMean),axis image, 
        end
        h = colorbar; title(h,'          10^')
        caxis([prctile(Var.StimMean(:),opt.ampLimit(1)) prctile(Var.StimMean(:),opt.ampLimit(2))])
        axis off
        title('Variance across stimuli \newline (averaged across reps)','fontsize',16)
        
        figure(fig2), subplot(1,3,2), 
        if opt.flip
            imagesc(fliplr(Var.StimRel)),axis image,
        else
            imagesc(Var.StimRel),axis image,  
        end
        h = colorbar; title(h,'          10^')
        caxis([prctile(Var.StimRel(:),opt.ampLimit(1)) prctile(Var.StimRel(:),opt.ampLimit(2))])
        axis off
        title('Normalized variance across stimuli \newline (averaged across reps)', 'fontsize',16)
    end
end
%% averaged variance of frames, across stim and reps
Var.FrameMean =   squeeze(mean(Var_frame,1));
% temp = abs(squeeze(mean(img_relamp,1))); temp(temp<0) = NaN;
Var.FrameRel =    log10(Var.FrameMean./Mean_frame);    % the relative variance (normalized by mean amplitude)
Var.FrameMean = log10(Var.FrameMean);

% Var.FrameMean =   log10(squeeze(mean(Var_frame,1)));
% Var.FrameRel =    log10(Var.FrameMean./abs(squeeze(mean(img_relamp,1))));    % the relative variance (normalized by mean amplitude)

if opt.plotON
    figure(fig1), subplot(1,3,3),
    if opt.flip
        imagesc(fliplr(Var.FrameMean)),axis image,
    else
        imagesc(Var.FrameMean),axis image, 
    end
    h = colorbar; title(h,'          10^')
    caxis([prctile(Var.FrameMean(:),opt.ampLimit(1)) prctile(Var.FrameMean(:),opt.ampLimit(2))])
    axis off
    title('Variance across frames \newline (averaged across reps and stimulus)', 'fontsize',16)
    
    figure(fig2), subplot(1,3,3),
    if opt.flip
        imagesc(fliplr(Var.FrameRel)),axis image,
    else
        imagesc(Var.FrameRel),axis image,  
    end
    h = colorbar; title(h,'          10^')
%     caxis([-2.5 -1.2])
    caxis([prctile(Var.FrameRel(:),opt.ampLimit(1)) prctile(Var.FrameRel(:),opt.ampLimit(2))])
    axis off
    title('Normalized variance across frames \newline (averaged across reps and stimulus)', 'fontsize',16)
end
end