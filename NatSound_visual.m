%% ==== visualize selected trials ====
% clear all
% load('\\10.16.58.229\Test_Imaging\2018.11.T1 (Marmoset 132D, Xintrinsic)\M132D-190119-10\190119T114434_Blue_Koehler_Fluo_GFP_P1.mat')
load('D:\Dropbox\_Yueqi Guo\_research\Imaging\=data=\data_wf\180724T140401_Blue_Koehler_Fluo_GFP_P1.mat')
% load('\\10.16.58.229\Test_Imaging\2018.03 T2 (Marmoset 80Z, Xintrinsic, Green & Fluo)\M80Z-180724-13\180724T140401_Blue_Koehler_Fluo_GFP_P1.mat')
%%
para.nRep = size(P.ProcDataMat,1);
para.nStim = size(P.ProcDataMat,2);
para.height = size(P.ProcDataMat,3);
para.width = size(P.ProcDataMat,4);
para.nFrame = size(P.ProcDataMat,5);

%P.ProcDataMat dimensions: (rep, trial, x(75), y(120), time(5*5))

para.fr = 5; 
para.preStim = 1;
para.durStim = 1;
para.postStim = 3;


%% plot a selected trial, all reps separately
iTrial = 31;
% ==== calculate relative amplitude change ====
mov = squeeze(P.ProcDataMat(:,iTrial,:,:,:)); 
mov_rel = mov;
img_rel = zeros(para.nRep+1, para.height, para.width); % mean image during response window 

for i = 1:para.nRep
    mov_temp = squeeze(mov(i,:,:,:)); % 75 x 120 x 25
    img_base = squeeze(mean(mov_temp(:,:,1:floor(para.fr*para.preStim)),3));
    img_base = repmat(img_base,1,1,para.nFrame);
    mov_rel(i,:,:,:) = (mov_temp - img_base)./img_base;
    img_rel(i,:,:) = mean(mov_rel(i,:,:,floor(para.fr*para.preStim)+1:end),4);
end

% ==== calculate mean movie ==== 
mov_mean = squeeze(mean(mov));
img_base = squeeze(mean(mov_mean(:,:,1:floor(para.fr*para.preStim)),3));
img_base = repmat(img_base,1,1,para.nFrame);
mov_rel(i+1,:,:,:) = (mov_mean - img_base)./img_base;
img_rel(i+1,:,:) = mean(mov_rel(i+1,:,:,floor(para.fr*para.preStim):end),4);


% ==== show movie as combined matrix (the last is the averaged response)
figure,
img_all = reshape(permute(img_rel,[2,3,1]),para.height,para.width*(para.nRep+1));

h = imagesc(img_all,[0,0.35]);
colorbar; 
axis image

pause
mov_all = zeros(para.height, para.width*(para.nRep+1));
for i = 1:para.nFrame
    mov_all = mov_rel(:,:,:,i);
    mov_all = reshape(permute(mov_all,[2,3,1]),para.height,para.width*(para.nRep+1));
    set(h,'CData',mov_all)
    pause(0.5)
end
set(h,'CData',img_all)

%% save as a video


%% show movies as separate matrices
figure;
for iRep = 1:para.nRep
    temp = squeeze(mov_rel(iRep,:,:,:));
    subplot(2,3,iRep); h(iRep) = imagesc(temp(:,:,1),[0,0.35]);colorbar;
    axis image
end
pause
for i = 1:para.nFrame
    for iRep = 1:para.nRep
    set(h(iRep),'CData',squeeze(mov_rel(iRep,:,:,i)))
    end
    pause(1)
end


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