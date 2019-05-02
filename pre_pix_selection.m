%% preprocessing - pixel selection
%% load pre-processed data
clear all
% load('\\10.16.58.229\Test_Imaging\2018.03 T2 (Marmoset 80Z, Xintrinsic, Green & Fluo)\M80Z-180724-13\180724T140401_Blue_Koehler_Fluo_GFP_P1.mat')
% load('\\10.16.58.229\Test_Imaging\2018.11.T1 (Marmoset 132D, Xintrinsic)\M132D-190119-10\190119T114434_Blue_Koehler_Fluo_GFP_P1.mat')
load('D:\Dropbox\_Yueqi Guo\_research\Imaging\=data=\data_wf\190119T114434_Blue_Koehler_Fluo_GFP_P1.mat')
% read background image
I = imread('D:\Dropbox\_Yueqi Guo\_research\Imaging\=data=\data_wf\190119T132947_Blue_Koehler_Fluo_GFP','tif');
I = double(imresize(I,1/16));
%
% Average across reps
nRep = size(P.ProcDataMat,1);
nTrial = size(P.ProcDataMat,2);
height = size(P.ProcDataMat,3);
width = size(P.ProcDataMat,4);
nFrame = size(P.ProcDataMat,5);
nStim = size(P.ProcDataMat,2);

%P.ProcDataMat dimensions: (rep, trial, x(75), y(120), time(5*5))

fr = 5; 
preStim = 1;
durStim = 2;
postStim = 2;


%% Calculate response variance across reps

V = zeros(nRep,nStim,width*height);
for i = 1:nRep
    mov = squeeze(P.ProcDataMat(i,:,:,:,:));

    img_base = squeeze(mean(mov(:,:,:,1:floor(fr*preStim)),4)); % first second: baseline
    img_amp = squeeze(mean(mov(:,:,:,floor(fr*preStim)+1:floor(fr*preStim+ + fr*durStim)),4));
    % figure,imshow(img_base,[])
    % figure,imshow(img_amp,[])
    img_relamp = (img_amp - img_base)./img_base;
    V(i,:,:) = reshape(img_relamp,nStim,width*height);
end

I_var = var(V);
I_var = reshape(squeeze(I_var),nStim,height,width);
%% see variance map for each stimulus
for i = 11:20
figure,imagesc(squeeze(I_var(i,:,:))),axis image
end
%% averaged variance across stimulus
I_var_mean = mean(I_var,1);
figure,imagesc(squeeze(I_var_mean)),axis image
%%
I_var_rel = squeeze(I_var_mean)./squeeze(mean(img_relamp,1));
figure, imagesc(I_var_rel,[0 0.1]),axis image

