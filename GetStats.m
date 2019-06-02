function [Mean, Var] = GetStats(para,DataMat,mask_final)

Mean =  zeros(para.nStim,1);
Var =   zeros(para.nStim,1);
for iTrial = 1:para.nStim
    mov = DataMat(:,iTrial,:,:,:);
    % average across reps
    mov_mean =      squeeze(mean(mov,1));
    % calculate deltaF/F (movie)
    img_base =      squeeze(mean(mov_mean(:,:,1:floor(para.fr*para.preStim)),3)); % pre-stimulus: baseline        
    img_base =      repmat(img_base,1,1,para.nFrame);
    mov_rel =       (mov_mean - img_base)./img_base;
    
    Var_temp =      var(mov_rel,0,3);
    % calculate deltaF/F (averaged image)
    img_rel =       squeeze(mean(mov_rel(:,:,floor(para.fr*para.preStim)+1 : floor(para.fr*(para.preStim + para.durStim))),3));
    img_rel =       img_rel.*double(mask_final);
    Mean(iTrial) =  mean(img_rel(img_rel~=0),'all');
    Var(iTrial) =   mean(Var_temp,'all');
end


end
