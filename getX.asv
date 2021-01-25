function [X, mov_rel, mov_rel_sep] = getX(DataMat, para, opt)
% mov_rel_sep: same size as DataMat, separately for each repetition

% set default parameters for opt
if ~isfield(opt, 'reps')
    opt.reps = 1:para.nRep;
end
if ~isfield(opt, 'tWindow')
    opt.tWindow = [para.preStim, para.preStim + para.durStim];
end
if ~isfield(opt, 'trials')
    opt.trials = 1:para.nStim;
end

 
nPanels = length(opt.trials);
mov_rel = zeros(nPanels,para.height,para.width,para.nFrame);
mov_rel_sep = zeros(length(opt.reps), nPanels,para.height,para.width,para.nFrame);
X = zeros(nPanels,para.width*para.height);

% reps = 1:para.nRep;
for i = 1:nPanels
    iTrial = opt.trials(i);
    mov = DataMat(opt.reps,iTrial,:,:,:);
    % for each repetition
    for j = 1:length(opt.reps)
        iRep = opt.reps(j);
        img_base1 = squeeze(mean(mov(j,:,:,:,1:floor(para.fr*para.preStim)),5)); %changed iRep to j within mov
        img_base1 = repmat(img_base1,1,1,para.nFrame);
        mov_rel_sep(j,i,:,:,:) = squeeze(mov(j,:,:,:,:))./img_base1;
    end
    
    mov_mean = squeeze(mean(mov,1)); % height x width x frames
    img_base = squeeze(mean(mov_mean(:,:,1:floor(para.fr*para.preStim)),3)); % first second: baseline
    img_base = repmat(img_base,1,1,para.nFrame);
    mov_rel(i,:,:,:) = (mov_mean - img_base)./img_base;
    % response window: same as stimulus period
%     img_amp = squeeze(mean(mov_mean(:,:,floor(para.fr*para.preStim)+1:floor(para.fr*para.preStim+ + para.fr*para.durStim + 8)),3));
    img_rel(i,:,:) = squeeze(mean(mov_rel(i,:,:,floor(para.fr*opt.tWindow(1))+1 : floor(para.fr*opt.tWindow(2))),4));  
end
X = reshape(img_rel, nPanels, para.width*para.height);
end