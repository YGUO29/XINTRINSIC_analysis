function X = getX(DataMat, para, opt)
if isempty(opt.reps)
    opt.reps = 1:para.nRep;
end
if isempty(opt.tWindow)
    opt.tWindow = [para.preStim, para.preStim + para.durStim];
end

X = zeros(para.nStim,para.width*para.height);
% reps = 1:para.nRep;
for i = 1:para.nStim
    iTrial = opt.trials(i);
    mov = DataMat(opt.reps,iTrial,:,:,:);
    mov_mean = squeeze(mean(mov,1)); % height x width x frames
    img_base = squeeze(mean(mov_mean(:,:,1:floor(para.fr*para.preStim)),3)); % first second: baseline
    
    img_base = repmat(img_base,1,1,para.nFrame);
    mov_rel(i,:,:,:) = (mov_mean - img_base)./img_base;
    % response window: same as stimulus period
%     img_amp = squeeze(mean(mov_mean(:,:,floor(para.fr*para.preStim)+1:floor(para.fr*para.preStim+ + para.fr*para.durStim + 8)),3));
    img_rel(i,:,:) = squeeze(mean(mov_rel(i,:,:,floor(para.fr*opt.tWindow(1))+1 : floor(para.fr*opt.tWindow(2))),4));  
end
X = reshape(img_rel, para.nStim, para.width*para.height);
end