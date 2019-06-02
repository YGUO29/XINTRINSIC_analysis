clear all

MonkeyID = 2; % 1: 80Z; 2: 132D session 1; 3: 132D session 2;

if MonkeyID == 1 % 80Z
file_mat = '180724T140401_Blue_Koehler_Fluo_GFP_P1.mat';
file_tif = '180724T133546_Blue_Koehler_Fluo_GFP.tif';
path_mat = '\\FANTASIA-DS3617\Test_Imaging\2018.03 T2 (Marmoset 80Z, Xintrinsic, Green & Fluo)\M80Z-180724-13\';
    load(fullfile(path_mat,file_mat));
    [~,nametemp,~,] = fileparts(file_mat);
    load(fullfile(path_mat,[nametemp(1:35),'.mat']));
    I = imread(fullfile(path_mat,file_tif));
    
elseif MonkeyID == 2 % 132D 1st natural sound
file_mat = '190119T114434_Blue_Koehler_Fluo_GFP_P1.mat';
file_tif = '190119T132947_Blue_Koehler_Fluo_GFP.tif';
path_mat = '\\FANTASIA-DS3617\Test_Imaging\2018.11.T1 (Marmoset 132D, Xintrinsic)\M132D-190119-10\';
    load(fullfile(path_mat,file_mat));
    [~,nametemp,~,] = fileparts(file_mat);
    load(fullfile(path_mat,[nametemp(1:end-3),'.mat']));
    I = imread(fullfile(path_mat,file_tif));
    
elseif MonkeyID == 3 % 132D 2nd natural sound
file_mat = '190224T095545_Blue_Koehler_Fluo_GFP_P1.mat';
file_tif = '190224T121936_Blue_Koehler_Fluo_GFP.tif';
path_mat = '\\FANTASIA-DS3617\Test_Imaging\2018.11.T1 (Marmoset 132D, Xintrinsic)\M132D-190224-09\';
    load(fullfile(path_mat,file_mat));
    [~,nametemp,~,] = fileparts(file_mat);
    load(fullfile(path_mat,[nametemp(1:end-3),'.mat']));
    I = imread(fullfile(path_mat,file_tif));

elseif MonkeyID == 4 % 132D 1st scrambled sound
load('\\FANTASIA-DS3617\Test_Imaging\2018.11.T1 (Marmoset 132D, Xintrinsic)\M132D-190409-10\190409T133916_Blue_Koehler_Fluo_GFP_P1.mat')
load('\\FANTASIA-DS3617\Test_Imaging\2018.11.T1 (Marmoset 132D, Xintrinsic)\M132D-190409-10\190409T133916_Blue_Koehler_Fluo_GFP.mat')
I = imread('\\FANTASIA-DS3617\Test_Imaging\2018.11.T1 (Marmoset 132D, Xintrinsic)\M132D-190409-10\190409T141027_Blue_Koehler_Fluo_GFP','tif');    
elseif MonkeyID == 5
load('\\FANTASIA-DS3617\Test_Imaging\2018.11.T1 (Marmoset 132D, Xintrinsic)\M132D-190419-09\190419T105449_Blue_Koehler_Fluo_GFP_P1.mat')
load('\\FANTASIA-DS3617\Test_Imaging\2018.11.T1 (Marmoset 132D, Xintrinsic)\M132D-190419-09\190419T105449_Blue_Koehler_Fluo_GFP.mat')
I = imread('\\FANTASIA-DS3617\Test_Imaging\2018.11.T1 (Marmoset 132D, Xintrinsic)\M132D-190419-09\190419T122732_Blue_Koehler_Fluo_GFP','tif');        
else 
        [file_mat,path_mat] = uigetfile('*.mat','Select a mat file to analyze','\\FANTASIA-DS3617\Test_Imaging\2018.11.T1 (Marmoset 132D, Xintrinsic)');
        load(fullfile(path_mat,file_mat));clc
        [~,nametemp,~,] = fileparts(file_mat);
        load(fullfile(path_mat,[nametemp(1:36),'.mat']));
        [file_tif,path_tif] = uigetfile('*.tif','Select a surface image','\\FANTASIA-DS3617\Test_Imaging\2018.11.T1 (Marmoset 132D, Xintrinsic)');
        I = imread(fullfile(path_tif,file_tif));
end

%%
I = double(imresize(I,1/16));
para.nRep =     size(P.ProcDataMat,1);
para.nStim =    size(P.ProcDataMat,2);
para.height =   size(P.ProcDataMat,3);
para.width =    size(P.ProcDataMat,4);
para.nFrame =   size(P.ProcDataMat,5);


para.fr =       P.ProcFrameRate; 
para.preStim =  S.TrlDurPreStim;
para.durStim =  S.TrlDurStim;
para.postStim = S.TrlDurPostStim;
para.filename = file_mat;
para.pathname = path_mat;
DataMat =       P.ProcDataMat;

%% select trials to view (temp)
folder_sound = 'D:\=code=\McdermottLab\sound_natural\';
list = dir(fullfile(folder_sound,'*.wav'));
names_sound = natsortfiles({list.name})';
folder_sound = '\\FANTASIA-DS3617\Test_Imaging\=Sounds=\Natural_XINTRINSIC';
list = dir(fullfile(folder_sound,'*.wav'));
names_sound2 = natsortfiles({list.name})';
opt.trials = [];
for i = 1:length(names_sound2)-1
    opt.trials(i) = find(strcmp(names_sound2(i+1),names_sound));
end

%% View data
opt.ampLimit =  [-0.1 0.1];
opt.mode =      'avgrep'; % avgrep or allrep
opt.plotMode =  'combined';
% opt.trials =    1:para.nStim;
opt.trials =    1:36;
opt.saveON =    0;

ViewData(para,DataMat,opt)
%%  Variance across reps 1
opt.plotON =    1;
opt.ampLimit =  0.1;
[Var] =         AnalysisVar(para,DataMat,opt);
%% Generate masks for pixel selection
DataMat = P.ProcDataMat;
p_th = 0.0005; r_th = 0.3; c_th = 0.3;
plot_on = 1; 
[mask_sig, mask_consis, mask_corr] = PixSelect(para,DataMat,I,p_th,r_th,c_th,plot_on);
% ==== generate final masks with surface image ====
mask_final = mask_corr.*mask_sig;

%% Manually drawn circular mask
figure, imshow(I,[])
h = images.roi.Circle(gca,'Center',[floor(para.width/2) floor(para.height/2)],'Radius',40); 
mask_final = createMask(h);

%% get response Mean & Variance within ROI
[R_mean, R_var] = GetStats(para,DataMat,mask_final);
R_mean = R_mean.*100; % unit: %

folder_sound =                      'D:\=code=\McdermottLab\sound_natural\';
list =                              dir(fullfile(folder_sound,'*.wav'));
SoundName =                         natsortfiles({list.name})';
[R_mean_inorder,I_mean] =           sort(R_mean,'descend');
[R_var_inorder,I_var] =             sort(R_var,'descend');
[R_relvar_inorder,I_relvar] =       sort(R_var./R_mean,'descend');
[~,Rank_mean] =                     sort(I_mean);
[~,Rank_var] =                      sort(I_var);
[~,Rank_relvar] =                   sort(I_relvar);

SoundName_inorder =                 SoundName(I_relvar);


%% ============ construct a matrix for ICA analysis ===================
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
tic,[R, W] = nonparametric_ica(-X, K, RANDOM_INITS, PLOT_FIGURES);toc

comp = cell(1,K);
I_norm = (I - min(min(I)))./(max(max(I)) - min(min(I)));

figure,
[p,n] = numSubplots(K);
cutoff = 0.1;
for i = 1:K
    comp{i} = reshape(W(i,:),para.height,para.width);
%     subplot(p(1),p(2),i),imagesc(comp{i},cutoff.*[-1 1]),axis image
    mask = comp{i}; mask(mask > cutoff) = cutoff; mask(mask < - cutoff) = - cutoff;
    mask = mask.*8; 
    img = repmat(I_norm,1,1,3); % three layers, representing R,G,B 
    img(:,:,1) = img(:,:,1) + mask;
    subplot(p(1),p(2),i),imagesc(img),axis image
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
R_mean = R_inorder(:,i); index = I_inorder(:,i);
tags_inorder = tags(index);
Color_inorder = Color(tags_inorder,:);
subplot(3,1,i)
b = bar(R_mean,'FaceColor','flat');
title(['Response Magnitude, Component ',num2str(i)],'fontsize',16)
b.CData = Color_inorder;
end
figure,
for i = 4:6
% ==== select a component to look at
% [resp,index] = sort(R(:,i),'descend');
R_mean = R_inorder(:,i); index = I_inorder(:,i);
% arrange tags according to response profile
tags_inorder = tags(index);
Color_inorder = Color(tags_inorder,:);
subplot(3,1,i-3)
b = bar(R_mean,'FaceColor','flat');
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