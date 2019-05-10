


%% plot selected trials, avraged across all reps
% close all
figure,
trials = 30.*ones(1,6) + [1 2 3 4 5 6]; 
for i = 1:6
    iTrial = trials(i);
    mov = P.ProcDataMat(:,iTrial,:,:,:);
    mov_mean = squeeze(mean(mov,1));

    img_base = squeeze(mean(mov_mean(:,:,1:floor(para.fr*para.preStim)),3)); % first second: baseline
    img_amp = squeeze(mean(mov_mean(:,:,floor(para.fr*para.preStim)+1:floor(para.fr*para.preStim+ + para.fr*para.durStim)),3)); 
    % figure,imshow(img_base,[])
    % figure,imshow(img_amp,[])
    img_relamp = (img_amp - img_base)./img_base;

    subplot(1,6,i)
    imagesc(img_relamp,[0,0.3])
    axis image
    colorbar
end

%%
iTrial = 1;
mov = P.ProcDataMat(:,iTrial,:,:,:);
mov_mean = squeeze(mean(mov,1));

img_base = squeeze(mean(mov_mean(:,:,1:floor(para.fr*para.preStim)),3)); % first second: baseline
img_base = repmat(img_base,1,1,para.nFrame);
mov_rel = (mov_mean - img_base)./img_base;
img_all = squeeze(mean(mov_rel(:,:,para.fr*para.preStim:para.fr*(para.preStim+para.durStim)),3));

figure,
h = imagesc(img_all,[0,0.3])
axis image
colorbar
pause
mov_all = zeros(para.height, para.width*(para.nRep+1))


for i = 1:para.nFrame
    mov_all = mov_rel(:,:,i);
    set(h,'CData',mov_all)
    pause(0.2)
end
set(h,'CData',img_all)
