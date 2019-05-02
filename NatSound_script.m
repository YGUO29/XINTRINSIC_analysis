clear all
% load('\\10.16.58.229\Test_Imaging\2018.11.T1 (Marmoset 132D, Xintrinsic)\M132D-190119-10\190119T114434_Blue_Koehler_Fluo_GFP_P1.mat')
% load('\\10.16.58.229\Test_Imaging\2018.03 T2 (Marmoset 80Z, Xintrinsic, Green & Fluo)\M80Z-180724-13\180724T140401_Blue_Koehler_Fluo_GFP_P1.mat')
MonkeyID = 5; % 1: 80Z; 2: 132D session 1; 3: 132D session 2;

if MonkeyID == 1 % 80Z
load('\\FANTASIA-DS3617\Test_Imaging\2018.03 T2 (Marmoset 80Z, Xintrinsic, Green & Fluo)\M80Z-180724-13\180724T140401_Blue_Koehler_Fluo_GFP_P1.mat')
I = imread('\\FANTASIA-DS3617\Test_Imaging\2018.03 T2 (Marmoset 80Z, Xintrinsic, Green & Fluo)\M80Z-180724-13\180724T133546_Blue_Koehler_Fluo_GFP','tif');
elseif MonkeyID == 2 % 132D 1st natural sound
load('\\FANTASIA-DS3617\Test_Imaging\2018.11.T1 (Marmoset 132D, Xintrinsic)\M132D-190119-10\190119T114434_Blue_Koehler_Fluo_GFP_P1.mat')
I = imread('\\FANTASIA-DS3617\Test_Imaging\2018.11.T1 (Marmoset 132D, Xintrinsic)\M132D-190119-10\190119T132947_Blue_Koehler_Fluo_GFP','tif');        
elseif MonkeyID == 3 % 132D 2nd natural sound
load('\\FANTASIA-DS3617\Test_Imaging\2018.11.T1 (Marmoset 132D, Xintrinsic)\M132D-190224-09\190224T095545_Blue_Koehler_Fluo_GFP_P1.mat')
I = imread('\\FANTASIA-DS3617\Test_Imaging\2018.11.T1 (Marmoset 132D, Xintrinsic)\M132D-190224-09\190224T121936_Blue_Koehler_Fluo_GFP','tif');    
elseif MonkeyID == 4 % 132D 1st scrambled sound
load('\\FANTASIA-DS3617\Test_Imaging\2018.11.T1 (Marmoset 132D, Xintrinsic)\M132D-190409-10\190409T133916_Blue_Koehler_Fluo_GFP_P1.mat')
I = imread('\\FANTASIA-DS3617\Test_Imaging\2018.11.T1 (Marmoset 132D, Xintrinsic)\M132D-190409-10\190409T141027_Blue_Koehler_Fluo_GFP','tif');    
elseif MonkeyID == 5
load('\\FANTASIA-DS3617\Test_Imaging\2018.11.T1 (Marmoset 132D, Xintrinsic)\M132D-190119-10\190119T112530_Blue_Koehler_Fluo_GFP_P1.mat')
I = imread('\\FANTASIA-DS3617\Test_Imaging\2018.11.T1 (Marmoset 132D, Xintrinsic)\M132D-190119-10\190119T132947_Blue_Koehler_Fluo_GFP','tif');        
end
I = double(imresize(I,1/16));

para.nRep = size(P.ProcDataMat,1);
para.nStim = size(P.ProcDataMat,2);
para.height = size(P.ProcDataMat,3);
para.width = size(P.ProcDataMat,4);
para.nFrame = size(P.ProcDataMat,5);

%P.ProcDataMat dimensions: (rep, trial, x(75), y(120), time(5*5))

para.fr = 5; 
para.preStim = 1;
para.durStim = 2;
para.postStim = 2;
%%  Variance across reps 1
DataMat = P.ProcDataMat;
[Var_rep_mean, Var_rep_rel, Var_stim_mean, Var_stim_rel] = AnalysisVar(para,DataMat);
%% Generate masks for pixel selection
DataMat = P.ProcDataMat;
p_th = 0.0005; r_th = 0.3; c_th = 0.3;
plot_on = 1; 
[mask_sig, mask_consis, mask_corr] = PixSelect(para,DataMat,I,p_th,r_th,c_th,plot_on);
% ==== generate final masks with surface image ====
mask_final = mask_corr.*mask_sig;


%% construct a matrix for ICA analysis
X = zeros(para.nStim,para.width*para.height);
for i = 1:para.nStim
    mov = P.ProcDataMat(:,i,:,:,:);
    mov_mean = squeeze(mean(mov,1));

    img_base = squeeze(mean(mov_mean(:,:,1:floor(para.fr*para.preStim)),3)); % first second: baseline
    img_amp = squeeze(mean(mov_mean(:,:,floor(para.fr*para.preStim)+1:floor(para.fr*para.preStim+ + para.fr*para.durStim)),3));
    img_relamp = (img_amp - img_base)./img_base;
%     img_relamp = img_relamp.*mask_corr;
    X(i,:) = reshape(img_relamp,para.width*para.height,1);
end

%% perform ICA analysis
addpath('D:\=code=\McdermottLab\toolbox_nonparametric-ICA')
K = 6;
RANDOM_INITS = 10;
PLOT_FIGURES = 0;
tic,[R, W] = nonparametric_ica(X, K, RANDOM_INITS, PLOT_FIGURES);toc

comp = cell(1,K);

I_norm = (I - min(min(I)))./(max(max(I)) - min(min(I)));
figure,
for i = 1:K
    comp{i} = reshape(W(i,:),para.height,para.width);
%     subplot(3,3,i),imagesc(comp{i},0.15.*[-1 1]),axis image
    mask = comp{i}; mask(mask > 0.15) = 0.15; mask(mask < -0.15) = -0.15;
    mask = mask.*8; 
    img = repmat(I_norm,1,1,3); % three layers, representing R,G,B 
    img(:,:,1) = img(:,:,1) + mask;
    subplot(3,3,i),imagesc(img),axis image
end

%% analyze response amplitude for components
plot_on = 1;
[I_inorder, R_inorder, tags_inorder, snames_inorder] = CompRespCat(R,plot_on);
%% calculate percentage variance explained by this number of components
V = zeros(2,para.nStim,para.width*para.height);
ind = [1,2]; % select two reps to compare
for i = 1:2
    mov = squeeze(P.ProcDataMat(ind(i),:,:,:,:));

    img_base = squeeze(mean(mov(:,:,:,1:floor(para.fr*para.preStim)),4)); % first second: baseline
    img_amp = squeeze(mean(mov(:,:,:,floor(para.fr*para.preStim)+1:floor(para.fr*para.preStim+ + para.fr*para.durStim)),4));
    % figure,imshow(img_base,[])
    % figure,imshow(img_amp,[])
    img_relamp = (img_amp - img_base)./img_base;
    V(i,:,:) = reshape(img_relamp,para.nStim,para.width*para.height);
end

V1 = squeeze(V(1,:,:)); 
V2 = squeeze(V(2,:,:)); 

rho = zeros(para.width*para.height,1);
for i = 1:para.height*para.width
    v1 = V1(:,i); v2 = V2(:,i);
    v_proj1 = R*pinv(R)*v1;
    v_proj2 = R*pinv(R)*v2;
    rho1 = corr(v_proj1, v2);
    rho2 = corr(v_proj2, v1);
    rho(i) = tanh(0.5*(atanh(rho1)+atanh(rho2)));
%     rho(i) = rho(i)/sqrt(corr(v1,v2)*corr(v_proj1,v_proj2));
%     rho(i) = abs(rho(i)).^2;
%     rho(i) = corr(v_proj1,v_proj2).^2;
end
median(rho)

%% Tonotopic map analysis
X1 = X(1:7,:);
X2 = X(8:end,:);
freqs = logspace(log10(220),log10(14080),7);
% calculate weighted preferred frequency
I_bf1 = zeros(width*height,1);
figure
for i = 1:width*height
    temp = X1(:,i); temp(temp<0) = 0;
    if var(temp./sum(temp))>0
%         I_bf1(i) = freqs*(temp./sum(temp));
        [~,ind] = max(temp);
        I_bf1(i) = freqs(ind);
        plot(temp),hold on
    end

    
end
I_bf1 = reshape(I_bf1,height,width);
figure,imagesc(I_bf1),axis image

% I_bf2 = zeros(width*height,1);
% for i = 1:width*height
%     temp = X2(:,i); temp(temp<0) = 0;
%     I_bf2(i) = freqs*(temp./sum(temp));
% %     [~,ind] = max(temp);
% %     I_bf2(i) = freqs(ind);
% end
% I_bf2 = reshape(I_bf2,height,width);
% figure,imagesc(I_bf2),axis image

%% look at response and stimulus

% filelist = dir('\\10.16.58.229\Test_Procedures\==Slides & Documents\Music\naturalsounds165\naturalsounds165')
folder_origin = 'D:\Dropbox\_Yueqi Guo\_research\Imaging\=code=\McdermottLab\sound_natural';
% list = dir(fullfile(folder_origin,'*.wav'));
% Names = natsortfiles({list.name})';
tags = xlsread([folder_origin,'\NatSound_label'],1);
nTags = max(tags);
Color = zeros(nTags,3);
% set colors
Color(1,:) = [19  78  150]./255;
Color(2,:) = [0   153 211]./255;
Color(3,:) = [0   104 78 ]./255;
Color(4,:) = [48  191 159]./255;
Color(5,:) = [124 66  150]./255;
Color(6,:) = [162 121 186]./255;
Color(7,:) = [202 51  32 ]./255;
Color(8,:) = [255 103 132]./255;
Color(9,:) = [119 119 119]./255;
Color(10,:) = [255 139 24]./255;

%% bar plot figure
[R_inorder,I_inorder] = sort(R,'descend');

close all
figure,
for i = 1:3
% ==== select a component to look at
% [resp,index] = sort(R(:,i),'descend');
% arrange tags according to response profile
resp = R_inorder(:,i); index = I_inorder(:,i);
tags_inorder = tags(index);
Color_inorder = Color(tags_inorder,:);
subplot(3,1,i)
b = bar(resp,'FaceColor','flat');
title(['Response Magnitude, Component ',num2str(i)],'fontsize',16)
b.CData = Color_inorder;
end
figure,
for i = 4:6
% ==== select a component to look at
% [resp,index] = sort(R(:,i),'descend');
resp = R_inorder(:,i); index = I_inorder(:,i);
% arrange tags according to response profile
tags_inorder = tags(index);
Color_inorder = Color(tags_inorder,:);
subplot(3,1,i-3)
b = bar(resp,'FaceColor','flat');
title(['Response Magnitude, Component ',num2str(i)],'fontsize',16)
b.CData = Color_inorder;
end

%% plot response amplitude according to catagories
RespGroupMean = zeros(6,nTags);
RespGroupStd  = RespGroupMean;
for t = 1:nTags
    for i = 1:6
    %     resp = R_inorder(:,i); index = I_inorder(:,i);
    index = find(tags == t);
    resp_temp = R(index,i);
    RespGroupMean(i,t) = mean(resp_temp);
    RespGroupStd(i,t) = std(resp_temp);
    end
end

figure,
for i = 1:6
subplot(2,3,i),errorbar(RespGroupMean(i,:),RespGroupStd(i,:))
end