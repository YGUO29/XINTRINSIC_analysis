function X = getX(DataMat, para, opt)

X = zeros(para.nStim,para.width*para.height);
reps = 1:para.nRep;
for i = 1:para.nStim
    mov = DataMat(reps,i,:,:,:);
    mov_mean = squeeze(mean(mov,1)); % height x width x frames

    img_base = squeeze(mean(mov_mean(:,:,1:floor(para.fr*para.preStim)),3)); % first second: baseline
    
    % response window: same as stimulus period
    img_amp = squeeze(mean(mov_mean(:,:,floor(para.fr*para.preStim)+1:floor(para.fr*para.preStim+ + para.fr*para.durStim + 8)),3));

    % response window: after stimulus (offset)
%     img_amp = squeeze(mean(mov_mean(:,:,floor(para.fr*para.preStim+ + para.fr*para.durStim)+1:end),3));

    img_relamp = (img_amp - img_base)./img_base;
%     img_relamp = img_relamp.*mask_final;
%     img_relamp = img_relamp.*mask_corr;
    X(i,:) = reshape(img_relamp,para.width*para.height,1);
end