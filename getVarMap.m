function Var = getVarMap(para,DataMat,opt)
% DataMat = [rep, trial, height, width, frams]
% opt.plotON, opt.ampLimit

% Construct a matrix V to compute variance across REP and STIM
% [rep, stim, pixels(w*h)]
if ~isfield(opt,'tWindow')
    opt.tWindow = [para.preStim para.durStim + para.preStim];
end
if ~isfield(opt,'plotON')
    opt.plotON = 1;
end


V = zeros(para.nRep,para.nStim,para.width*para.height);
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

if para.nRep ~= 1
Var_rep = var(V);
Var_rep = reshape(squeeze(Var_rep),para.nStim,para.height,para.width); % [stim, height, width]
end

if para.nStim ~= 1
Var_stim = var(V,0,2);
Var_stim = reshape(squeeze(Var_stim),para.nRep,para.height,para.width); % [rep, height, width]
end

%% Compute variance across frames
Var_frame = zeros(para.nStim,para.height,para.width);
Mean_frame = zeros(para.nStim,para.height,para.width);
% Base_frame = zeros(para.nStim,para.height,para.width); % temporary

for i = 1:para.nStim
% for i = 1:10
    mov_temp = mean(squeeze(DataMat(:,i,:,:,:)),1); % 1 x height x width x frames
    mov_temp = squeeze(mov_temp); % height x width x frames
    img_base = squeeze(mean(mov_temp(:,:,1:floor(para.fr*para.preStim)),3));
%     Base_frame(i,:,:) = img_base; % temporary - keep base frames for each stimulus
    
    img_base = repmat(img_base,1,1,para.nFrame);
    mov_rel = (mov_temp - img_base)./img_base;
    Var_frame(i,:,:) = var(mov_rel(:,:,para.fr*opt.tWindow(1)+1:floor(para.fr*opt.tWindow(2))),0,3);
    Mean_frame(i,:,:) = mean(mov_rel(:,:,para.fr*opt.tWindow(1)+1:floor(para.fr*opt.tWindow(2))),3);
end

%% see variance map for each stimulus
% for i = 1:10
% figure,imagesc(squeeze(Var_frame(i,:,:))),axis image
% end
fig1 = figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w','position',[943 600 2304 420]);
fig2 = figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w','position',[943 600 2304 420]);

%% averaged variance of reps across reps
if para.nRep ~= 1
    Var.RepMean =   squeeze(mean(Var_rep,1));
    temp = abs(squeeze(mean(img_relamp,1))); temp(temp<0) = NaN;
    Var.RepRel =    log10(Var.RepMean./temp);    % the relative variance (normalized by mean amplitude)
    Var.RepMean = log10(Var.RepMean);

    if opt.plotON
        figure(fig1), subplot(1,3,1)
        imagesc(Var.RepMean),axis image, 
        h = colorbar; title(h,'          10^')
        caxis([prctile(Var.RepMean(:),opt.ampLimit(1)) prctile(Var.RepMean(:),opt.ampLimit(2))])
        axis off
        title('Variance across repetitions \newline (averaged across stimulus)','fontsize',16)
        
        figure(fig2), subplot(1,3,1)
        imagesc(Var.RepRel),axis image, 
        h = colorbar; title(h,'          10^')
        caxis([prctile(Var.RepRel(:),opt.ampLimit(1)) prctile(Var.RepRel(:),opt.ampLimit(2))])
        axis off
        title('Normalized variance across repetitionss \newline (averaged across stimulus)', 'fontsize',16)
    end
end
%% averaged variance of stims across stimuli
if para.nStim ~= 1
    Var.StimMean =   squeeze(mean(Var_stim,1));
    temp = abs(squeeze(mean(img_relamp,1))); temp(temp<0) = NaN;
    Var.StimRel =    log10(Var.StimMean./temp);    % the relative variance (normalized by mean amplitude)
    Var.StimMean = log10(Var.StimMean);
    
    if opt.plotON
        figure(fig1), subplot(1,3,2), 
        imagesc(Var.StimMean),axis image,  
        h = colorbar; title(h,'          10^')
        caxis([prctile(Var.StimMean(:),opt.ampLimit(1)) prctile(Var.StimMean(:),opt.ampLimit(2))])
        axis off
        title('Variance across stimuli \newline (averaged across reps)','fontsize',16)
        
        figure(fig2), subplot(1,3,2), 
        imagesc(Var.StimRel),axis image,  
        h = colorbar; title(h,'          10^')
        caxis([prctile(Var.StimRel(:),opt.ampLimit(1)) prctile(Var.StimRel(:),opt.ampLimit(2))])
        axis off
        title('Normalized variance across stimuli \newline (averaged across reps)', 'fontsize',16)
    end
end
%% averaged variance of frames, across stim and reps
 Var.FrameMean =   squeeze(mean(Var_frame,1));
temp = abs(squeeze(mean(img_relamp,1))); temp(temp<0) = NaN;
Var.FrameRel =    log10(Var.FrameMean./temp);    % the relative variance (normalized by mean amplitude)
Var.FrameMean = log10(Var.FrameMean);

% Var.FrameMean =   log10(squeeze(mean(Var_frame,1)));
% Var.FrameRel =    log10(Var.FrameMean./abs(squeeze(mean(img_relamp,1))));    % the relative variance (normalized by mean amplitude)

if opt.plotON
    figure(fig1), subplot(1,3,3),
    imagesc(Var.FrameMean),axis image,  
    h = colorbar; title(h,'          10^')
    caxis([prctile(Var.FrameMean(:),opt.ampLimit(1)) prctile(Var.FrameMean(:),opt.ampLimit(2))])
    axis off
    title('Variance across frames \newline (averaged across reps and stimulus)', 'fontsize',16)
    figure(fig2), subplot(1,3,3),
    imagesc(Var.FrameRel),axis image,  
    h = colorbar; title(h,'          10^')
    caxis([prctile(Var.FrameRel(:),opt.ampLimit(1)) prctile(Var.FrameRel(:),opt.ampLimit(2))])
    axis off
    title('Normalized variance across frames \newline (averaged across reps and stimulus)', 'fontsize',16)
end
end