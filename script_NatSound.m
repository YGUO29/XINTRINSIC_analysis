% clear all

MonkeyID = 3; % 1: 80Z; 2: 132D session11; 3: 132D session2; 4: 132D 1st Scrambled; 6: 132D 2nd Scrambled

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
elseif MonkeyID == 5 % 132D 2nd scrambled sound
load('\\FANTASIA-DS3617\Test_Imaging\2018.11.T1 (Marmoset 132D, Xintrinsic)\M132D-190419-09\190419T105449_Blue_Koehler_Fluo_GFP_P1.mat')
load('\\FANTASIA-DS3617\Test_Imaging\2018.11.T1 (Marmoset 132D, Xintrinsic)\M132D-190419-09\190419T105449_Blue_Koehler_Fluo_GFP.mat')
I = imread('\\FANTASIA-DS3617\Test_Imaging\2018.11.T1 (Marmoset 132D, Xintrinsic)\M132D-190419-09\190419T122732_Blue_Koehler_Fluo_GFP','tif');        
else % browse to open a file
        [file_mat,path_mat] = uigetfile('*.mat','Select a mat file to analyze','\\FANTASIA-DS3617\Test_Imaging\2018.11.T1 (Marmoset 132D, Xintrinsic)');
        load(fullfile(path_mat,file_mat));clc
        [~,nametemp,~,] = fileparts(file_mat);
        load(fullfile(path_mat,[nametemp(1:35),'.mat']));
        [file_tif,path_tif] = uigetfile('*.tif','Select a surface image','\\FANTASIA-DS3617\Test_Imaging\2018.11.T1 (Marmoset 132D, Xintrinsic)');
        I = imread(fullfile(path_tif,file_tif));
end

%
% I = double(imresize(I,1/16));
para.nRep =     size(P.ProcDataMat,1);
para.nStim =    size(P.ProcDataMat,2);
para.height =   size(P.ProcDataMat,3);
para.width =    size(P.ProcDataMat,4);
para.nFrame =   size(P.ProcDataMat,5);


para.fr =       P.ProcFrameRate; % frame rate
para.preStim =  S.TrlDurPreStim;
para.durStim =  S.TrlDurStim;
para.postStim = S.TrlDurPostStim;
para.filename = file_mat;
para.pathname = path_mat;
DataMat =       P.ProcDataMat;
I = double(imresize(I,[para.height,para.width]));
clear S P
%% select trials to compare Calcium and intrinsic imaging (temp)
% the entire set of natural sounds
folder_sound = 'D:\=code=\McdermottLab\sound_natural\';
list = dir(fullfile(folder_sound,'*.wav'));
names_sound = natsortfiles({list.name})';
% the subset of natural sounds for intrinsic imaging
folder_sound = '\\FANTASIA-DS3617\Test_Imaging\=Sounds=\Natural_XINTRINSIC2';
list = dir(fullfile(folder_sound,'*.wav'));
names_sound2 = natsortfiles({list.name})';
% pick out trials for the selected natural sounds
opt.trials = [];
for i = 1:length(names_sound2)-1
    opt.trials(i) = find(strcmp(names_sound2(i+1),names_sound));
end

%% View data
opt.ampLimit    =  [0 0.3];
opt.mode        =  'avgrep'; % avgrep or allrep
opt.plotMode    =  'separate'; % combined or separate (video saving is only available for 'combined' mode)
% opt.trials =    44+[1:12,23:28]; 
opt.trials      = 165-11:165;
opt.saveON      = 0; 
opt.soundON     = 0;
opt.reps        = [];
opt.tWindow     = []; % start and end of integration window for calculating response amplitude

ViewData(para,DataMat,opt)
%%  Variance across reps 1
opt.plotON      = 1;
opt.ampLimit    = 0.1;
opt.tWindow     = []; % select time window for analysis (default is the stimulus-on period)
Var           = AnalysisVar(para,DataMat,opt);
%% ============== Generate a mask ==============
%% Manually drawn circular mask
figure, imshow(I,[])
h = images.roi.Circle(gca,'Center',[floor(para.width/2) floor(para.height/2)],'Radius',floor(para.height/2)); 
% h = images.roi.Rectangle(gca,'Position',[floor(para.width/2)-floor(para.height/2), 1,...
%                                          floor(para.height), floor(para.height)]); 
% h = images.roi.Ellipse(gca,'Center',[floor(para.width/2) floor(para.height/2)],'Semiaxes',[40 20]); 
title('Press Enter after the position is adjusted')
pause
mask_manual = createMask(h);
mask_final = mask_manual;

%% masks based on response properties
p_th = 0.0001; % significance threshold
r_th = 0.3; % reliability threshold
c_th = 0.3; % correlation (across reps) threshold
plot_on = 1; % plot: 0 = no plots; 1 = plot masks; 2 = plot all figures
[mask_sig, mask_consis, mask_corr] = PixSelect(para,DataMat,I,p_th,r_th,c_th,plot_on);

% ==== generate final masks with surface image ====
mask_5d = repmat(mask_final,1, 1, para.nRep, para.nStim, para.nFrame);
mask_5d = permute(mask_5d,[3,4,1,2,5]);
tic,DataMat_masked = DataMat.mask_5d;toc
%% get response Mean & Variance within ROI
[R_mean, R_var] = getStats(para,DataMat,mask_final);
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

SoundName_inorder =                 SoundName(I_mean);


%% ============ construct a matrix for ICA analysis ===================
X = zeros(para.nStim,para.width*para.height);
reps = 1:para.nRep;
for i = 1:para.nStim
    mov = DataMat(reps,i,:,:,:);
    mov_mean = squeeze(mean(mov,1)); % height x width x frames

    img_base = squeeze(mean(mov_mean(:,:,1:floor(para.fr*para.preStim)),3)); % first second: baseline
    
    % response window: same as stimulus period
    img_amp = squeeze(mean(mov_mean(:,:,floor(para.fr*para.preStim)+1:floor(para.fr*para.preStim+ + para.fr*para.durStim)),3));

    % response window: after stimulus (offset)
%     img_amp = squeeze(mean(mov_mean(:,:,floor(para.fr*para.preStim+ + para.fr*para.durStim)+1:end),3));

    img_relamp = (img_amp - img_base)./img_base;
%     img_relamp = img_relamp.*mask_final;
%     img_relamp = img_relamp.*mask_corr;
    X(i,:) = reshape(img_relamp,para.width*para.height,1);
end

% perform ICA analysis
addpath('D:\=code=\McdermottLab\toolbox_nonparametric-ICA')
K = 6;

RANDOM_INITS = 10;
PLOT_FIGURES = 0;
tic,[R, W] = nonparametric_ica(X, K, RANDOM_INITS, PLOT_FIGURES);toc

comp = cell(1,K);
I_norm = (I - min(min(I)))./(max(max(I)) - min(min(I)));

figure,
[p,n] = numSubplots(K);
cutoff = 0.1;
for i = 1:K
    comp{i} = reshape(W(i,:),para.height,para.width);
%     subplot(p(1),p(2),i),imagesc(comp{i},cutoff.*[-1 1]),axis image, colorbar
    mask = comp{i}; mask(mask > cutoff) = cutoff; mask(mask < - cutoff) = - cutoff;
    mask = mask.*8; 
    img = repmat(I_norm,1,1,3); % three layers, representing R,G,B 
    img(:,:,1) = img(:,:,1) + mask;
    subplot(p(1),p(2),i),imagesc(img),axis image
end

%% analyze response amplitude for components
plot_on = 1;
[I_inorder, R_inorder, tags_inorder, snames_inorder] = getResponseProfile(R,plot_on);
%% regression with sound features
load('D:\=code=\Sound_analysis\F_yg_marm_full.mat') 
nFeat       = size(F.F_mat,1);
result_p = zeros(nFeat,K);
result_r = zeros(nFeat,K);
% ============= Regression ===============
for iComp = 1:K
    yy = R(:,iComp);
    for iFeat = 1:nFeat
        xx = F.F_mat(iFeat,:)';
%         [p,rsq,~] = RSquared(xx,yy);
%         result_p(iFeat,iComp) = p(1);
%         result_r(iFeat,iComp) = sqrt(rsq);
        rr = corrcoef(xx,yy);
        result_r(iFeat,iComp) = rr(1,2);
    end
end

% ===========================================
% plot feature correlations
scr_size = get(0,'ScreenSize');
scr_size = scr_size(3:4); %width and height
% f1 = figure,set(gcf,'position',[1,scr_size(2)./2,scr_size(1),scr_size(2)./3]);
% f2 = figure,set(gcf,'position',[1,1,scr_size(1),scr_size(2)./3]);
f = figure; set(gcf,'position',[1,1,scr_size]);
ind = [1 2 5 6 4 3];
% ind = [1 2 3 4 5 6]; % plot order
for iComp = 1:K
    % correlation with frequency power
%     figure(f1)    
    subplot(2,K,iComp),
    plot(F.FreqBounds(1:end-1), result_r(1:F.nFreq,ind(iComp)), 'linewidth',2,'Marker','x')
    hold on, plot([F.FreqBounds(1:2), F.FreqBounds(end-1)], [0 0 0],'linestyle','--','color','k')    
    set(gca,'xscale','log');
    ymax = max(max(  abs( result_r(1:F.nFreq,1:K) )  )); ylim([-ymax,ymax]);
    set(gca,'xtick',F.FreqBounds(1:end-1))
    set(gca,'xticklabels',arrayfun(@num2str,F.FreqBounds(1:end-1),'UniformOutput',false), 'fontsize',12)
    xlabel('Frequency','fontsize',14),
    title(['Correlation coefficient, component #',num2str(ind(iComp))],'fontsize',14)

%     figure(f2)
%     subplot(2,K,iComp + K), 
%     spectemp_r = reshape(result_r(F.nFreq+1:F.nFreq+F.nSpectemp,ind(iComp)), ...
%                         size(F.spectemp_mod,1), size(F.spectemp_mod,2));
%     cmax = max(max(  abs( result_r(F.nFreq+1 : F.nFreq+F.nSpectemp,1:K) )  )); 
%     imagesc(flipud(spectemp_r),[-cmax, cmax]),colorbar, colormap('jet')
%     set(gca,'yticklabels',arrayfun(@num2str,fliplr(F.spec_mod_rates),'UniformOutput',false), 'fontsize',12)
%     set(gca,'xticklabels',arrayfun(@num2str,F.temp_mod_rates,'UniformOutput',false), 'fontsize',12)
%     xlabel('Temporal modulation rate (cycles/s)','fontsize',14),
%     ylabel('Spectral modulation rate (cycles/octave)','fontsize',14),
%     title(['Correlation coefficient, component #',num2str(ind(iComp))],'fontsize',14)

%   correlation with the full spectrotemporal modulation power map
    subplot(2,K,iComp + K), 
    spectemp_r = reshape(result_r(F.nFreq+1:F.nFreq+F.nSpectemp_full,ind(iComp)), ...
                        size(F.spectemp_mod_full,1), size(F.spectemp_mod_full,2));
    cmax = max(max(  abs( result_r(F.nFreq+1 : F.nFreq+F.nSpectemp_full,1:K) )  )); 
    imagesc(flipud(spectemp_r),[-cmax, cmax]),colorbar, colormap('jet')
    set(gca,'yticklabels',arrayfun(@num2str,fliplr(F.spec_mod_rates),'UniformOutput',false), 'fontsize',12)
    set(gca,'xticklabels',arrayfun(@num2str,F.temp_mod_rates,'UniformOutput',false), 'fontsize',12)
    xlabel('Temporal modulation rate (cycles/s)','fontsize',14),
    ylabel('Spectral modulation rate (cycles/octave)','fontsize',14),
    title(['Correlation coefficient, component #',num2str(ind(iComp))],'fontsize',14)

%     saveas(f,['80Z_session1_component_',num2str(iComp),'.png'])

end
% saveas(f1,'132D_session1_FreqPower.png')
% saveas(f2,'132D_session1_SpecTemp.png')
saveas(f,'132D_session1_reg_reordered_full.png')


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
rho_norm = rho;
for i = 1:para.height*para.width
    v1 = V1(:,i); v2 = V2(:,i);
    v_proj1 = R*pinv(R)*v1;
    v_proj2 = R*pinv(R)*v2;
    rho1 = corr(v_proj1, v2);
    rho2 = corr(v_proj2, v1);
    rho(i) = tanh(0.5*(atanh(rho1)+atanh(rho2)));
    
    rho_norm(i) = rho(i)/sqrt(corr(v1,v2)*corr(v_proj1,v_proj2));
    rho_norm(i) = abs(rho_norm(i)).^2;
%     rho(i) = corr(v_proj1,v_proj2).^2;
end
result(1) = median(rho);
result(2) = median(rho_norm)

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


