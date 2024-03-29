% =============== load files to analyze ================
% analysis program for trial-based experiments, run this section to
% initialize
addpath(genpath(cd))
clear all
para.nRep = 0;  
% 80Z Calcium: file_mat = '180724T140401_Blue_Koehler_Fluo_GFP_P1.mat';
% 132D Calcium 1: file_mat = '190119T114434_Blue_Koehler_Fluo_GFP_P1.mat';
% 132D Calcium 2: file_mat = '190224T095545_Blue_Koehler_Fluo_GFP_P1.mat';
%% add more repetitions
% 1. run this section to load any number of files from the same folder
% 2. run this section again to select a different set of files... continue
% to stack different repetitions as in "DataMat"

[file_mat,path_mat] = uigetfile('*.mat','Select a mat file to analyze','U:\', 'Multiselect', 'on');
if contains(file_mat, 'GFP')
    modality = 'GFP';
elseif contains(file_mat, 'PBS')
    modality = 'PBS';
else
    disp('check if filename contains GFP or PBS')
end

if ~iscell(file_mat) % only load one file
        load(fullfile(path_mat,file_mat));
        file_parts = strsplit(file_mat,{'_','.'});
%         ind = find(strcmp(file_parts,'GFP')); % GFP or PBS ? 
        ind = find(strcmp(file_parts, modality));
        nametemp = strjoin(file_parts(1:ind),'_');
        load(fullfile(path_mat,[nametemp,'.mat']));
        
        if para.nRep == 0
            para.nRep =     size(P.ProcDataMat,1);
            para.nStim =    size(P.ProcDataMat,2);
            para.height =   size(P.ProcDataMat,3);
            para.width =    size(P.ProcDataMat,4);
            para.nFrame =   size(P.ProcDataMat,5);

            para.fr =       P.ProcFrameRate; % frame rate
            para.preStim =  S.TrlDurPreStim;
            para.durStim =  S.TrlDurStim;
            para.postStim = S.TrlDurPostStim;
            para.order =    S.SesTrlOrderVec;
            DataMat =       P.ProcDataMat; % DataMat = [rep, trial, height, width, frams]
            para.filename = {file_mat}; % make it a cell array to be consistent in data structure
            para.pathname = {path_mat};
        else
%             DataMat(para.nRep+1 : para.nRep+size(P.ProcDataMat,1),:,:,:,:) = squeeze(P.ProcDataMat(:,:,6:95,:,:));
            DataMat(para.nRep+1 : para.nRep+size(P.ProcDataMat,1),:,:,:,:) = squeeze(P.ProcDataMat);
            para.nRep = para.nRep + size(P.ProcDataMat, 1);
        end
%         clear S P
else % load multiple files
    for i = 1:length(file_mat)
        load(fullfile(path_mat,file_mat{i}));
        file_parts = strsplit(file_mat{i},{'_','.'});
%         ind = find(strcmp(file_parts,'P1'));
        ind = find(strcmp(file_parts, modality));
        nametemp = strjoin(file_parts(1:ind),'_');
        load(fullfile(path_mat,[nametemp,'.mat']));

        if para.nRep + i == 1 % read the initial file, initialize parameters
        %     para.nRep =     size(P.ProcDataMat,1);
            para.nStim =    size(P.ProcDataMat,2);
            para.height =   size(P.ProcDataMat,3);
            para.width =    size(P.ProcDataMat,4);
            para.nFrame =   size(P.ProcDataMat,5);
            para.fr =       P.ProcFrameRate; % frame rate

            para.preStim =  S.TrlDurPreStim;
            para.durStim =  S.TrlDurStim;
            para.postStim = S.TrlDurPostStim;
            para.order(para.nRep+i,:) =    S.SesTrlOrderVec;
            DataMat =       P.ProcDataMat; % DataMat = [rep, trial, height, width, frams]
            para.filename = file_mat;
            para.pathname = path_mat;
        else % continuous reading 
            DataMat(para.nRep+1 : para.nRep+size(P.ProcDataMat,1),:,:,:,:) = squeeze(P.ProcDataMat);
        end
        
%     clear S P
    para.nRep = size(DataMat,1);
    
    disp([num2str(para.nRep),' finished'])
    end     
end
para.side = S.MkySide;
size(DataMat)
para
%% read saved sessions
% [file_mat,path_mat] = uigetfile('*.mat','Select a mat file to analyze','U:\', 'Multiselect', 'on');
clear opt_read
opt_read.animal = '102D'; 
opt_read.session = 'NatVocMM';
opt_read.modal = 'Calcium';
opt_read.date = '201226';
% opt_read.subset = 'Nat'; % only read subset of the data (comment this line if reading all data)
% 80Z tone session: '50dB'
[DataMat, para] = readSession(opt_read);
para.mirror = 1;

DataMat = double(DataMat);
opt = struct;
% opt.tWindow = [para.preStim.*ones(para.nStim,1), para.preStim.*ones(para.nStim,1) + dur_mat'];
% opt.tWindow = [para.preStim.*ones(para.nStim,1), para.preStim.*ones(para.nStim,1) + repelem(dur_mat,2)'];

[X, DataMat_norm, mov_rel_sep] = getX(DataMat, para, opt);
% generate X matrix for each repetition (used in noise correction)
X_sep = squeeze(mean(mov_rel_sep(:,:,:,:,floor(para.preStim*para.fr) +1 : floor((para.preStim+para.durStim)*para.fr)),5));
X_sep = reshape(X_sep, para.nRep, para.nStim, para.height*para.width); 

% ======== construct a X without NaN (X may contain NaNs if there are masked pixels) ========
[~, ind_delete]     = find( isnan(X) ); % linear index
para.ind_save       = setdiff(1:para.width*para.height, ind_delete);
X(:, ind_delete)    = [];

eval(['para_', opt_read.session, '_', opt_read.animal, '= para;'])
eval(['X_', opt_read.session, '_', opt_read.animal, '= X;'])
eval(['X_sep_', opt_read.session, '_', opt_read.animal, '= X_sep;'])

%% get tonotopy structure
% opt_read.animal = '80Z';
eval(['para = para_', opt_read.session, '_', opt_read.animal,';'])
eval(['X = X_', opt_read.session, '_', opt_read.animal,';'])


para = getTonotopyParameters(para, opt_read.animal);
T = getTonotopy(X, para);

% register if necessary
T = getRegisteredTonotopy(T, para, opt_read.animal, 'manual'); % for 132D: use manual alignment

% get tonotopy contour
M.saturation = T.sat_map; % used for masking
M.hue = T.hue_map; % used for plotting contour
M.map_rgb = T.tonotopy{3}; % chose a map to overlay the contour on
para.sat_lim = 0.02;
para.zindex = T.zindex;
ct = getContour(M, para);

% verification - plot tonotopy with contour
figurex; 
imagesc(M.map_rgb); 
% imagesc(map_hue); 
axis image, axis off
plotContour(ct, 'gray')
% set(gca,'XDir','reverse')

eval(['T_', opt_read.animal, '= T;'])

%% View data
DataMat = double(DataMat);
opt = struct;

% ==== simplified, just get data matrix X and DataMat_norm ====
% opt.trials = 162; 
% opt.reps = ;
% opt.tWindow     = [para.preStim, para.preStim + para.durStim + 4]; % start and end of integration window for calculating response amplitude
% [X, DataMat_norm, mov_rel_sep] = getX(DataMat, para, opt);

% ==== get data matrix X and view/save the video ====
opt.ampLimit    = 0.2*[-1 1];
% opt.reps          = [1:11+3, 11+6:para.nRep];
% opt.trials        = 1:36; 
% % intrinsic natural sound: integrate until 4s after sound offset
% opt.tWindow       = [para.preStim, para.preStim + para.durStim + 4]; 
% % vocalization sessions: variable trial durations
% opt.tWindow       = [para.preStim.*ones(para.nStim,1), para.preStim.*ones(para.nStim,1) + dur_mat'];

% % delete some trials/reps if it has too many movements
% opt.omit_trials   = [12 16 30 36 48 48 16 16]; 
% opt.omit_reps     = [1 1 1 1 1 5 7 10];

% opt.p             = [1 9]; % subplot rows and columes; 
% opt.mode          = 'allrep';
% opt.plotMode      = 'separate';
% opt.saveON        = 1; % if see error when saving video (with videoFileWriter), check if the video file name is a string 
% opt.soundON       = 1;
% opt.sounddata     = S.SesSoundWave;
% opt.sounddata = data;
% opt.color         = {'div', 'PRGn'}; % jet or {'div', 'PRGn'} or 'green' (no negative values)
% para.pathname     = 'D:\SynologyDrive\=data=\Video\';
% para.sessionname  = ['102D_210108_music_RiB_cutoff=', num2str(opt.ampLimit(2))];
    
f = figurex;
tic
[X, DataMat_norm, ~] = ViewData(DataMat, para, opt); 
toc

% ======== construct a X without NaN (X may contain NaNs if there are masked pixels) ========
[~, ind_delete]     = find( isnan(X) ); % linear index
para.ind_save       = setdiff(1:para.width*para.height, ind_delete);
X(:, ind_delete)    = [];


%% save MATLAB files
% [X, DataMat_norm, mov_rel_sep] = getX(DataMat, para, opt);
% DataMat = DataMat(:,1:2:end,:,:,:);
% DataMat_norm = DataMat_norm(1:2:end,:,:,:);
% X = X(1:2:end,:);
% para.nStim = 180;

animal = '102D'; 
session = 'NatVocMM';
modal = 'Calcium';
date = '201226';
datapath = 'E:\Dropbox\_Yueqi Guo\_research\Imaging\=data=\';
% save([datapath, '\', animal, '\DataMatReg_', animal, '_', session, '_', num2str(para.nRep), 'reps.mat'],...
%     'DataMat_reg', 'para', '-v7.3')
% save([datapath, '\', animal, '\DataMat_', modal, '_', date, '_', animal, '_', session, '_', num2str(para.nRep), 'reps_'...
%     num2str(para.height), 'x', num2str(para.width), 'x', num2str(para.fr), 'fps.mat'],...
%     'DataMat', 'para', '-v7.3')
% save([datapath, '\', animal, '\DataMatNorm_', modal, '_', date, '_', animal, '_', session, '_', num2str(para.nRep), 'reps_'...
%     num2str(para.height), 'x', num2str(para.width), 'x', num2str(para.fr), 'fps.mat'],...
%     'DataMat_norm', 'para', '-v7.3')
save([datapath, '\', animal, '\DataMatProc_', modal, '_', date, '_', animal, '_', session, '_', num2str(para.nRep), 'reps_'...
    num2str(para.height), 'x', num2str(para.width), 'x', num2str(para.fr), 'fps.mat'],...
    'X', 'para', '-v7.3')

% ======= save as python format =======
% data = permute(DataMat_norm,[2,3,4,1]);
% save(['python_', file_mat], 'data')
%% Get a response map 
% (for all stimuli in this set)
opt_resp         = struct;
opt_resp.method  = 'var';
opt_resp.title   = 'Activation map for ripples';
opt_resp.X       = X;
opt_resp.X(opt_resp.X <= 0) = nan;
opt_resp.ampLimit = [-6 -2]; % for var
% opt_resp.ampLimit = [-2.8 -0.8]; % for mean


X_tone = getRespMap(DataMat_norm, para, opt_resp);
set(gca,'XDir','reverse')
%%  Variance across reps/stims/frames
opt = struct;
opt.ampLimit    = [0 95]; % the percentile of the data distribution as upper limit
% 98 for 102D, 80Z; 95 for 132D
% opt.tWindow     = [4 5]; % seconds
opt.flip = 0;
Var = getVarMap(para, DataMat, opt);
%% Run decomposition
eval(['para = para_', opt_read.session, '_', opt_read.animal,';'])
eval(['X = X_', opt_read.session, '_', opt_read.animal,';'])

opt_decomp.fluo = 1; 
opt_decomp.method = 'mICA'; % 'mICA' or 'NMF' or 'PCA'
opt_decomp.nRows = 1;
opt_decomp.plotON = 1;
opt_decomp.Ks = 6;
para.mirror = 1;
% para.ct = ct(1);

% opt.X_test = X_test; % to calculate reconstruction error with
figurex;
Decomp_102D = runDecomp(X, opt_decomp, para);

% plotComp(1, Decomp_102D, opt_decomp, para)
plotTonotopyWithComponents(1, Decomp_102D, T_102D, para)

if strcmp(opt_decomp.method, 'PCA')
    % scree plot
    figurex; plot(cumsum(recon_error)./sum(recon_error))
    ylabel('Variance explained')
    xlabel('Number of components')
end

%% Match components among subjects  (use 102D K=6 as a reference)
% match based on response profiles
% concatenate all components
matched_comp = runMatchSubjects(Decomp_102D, Decomp_132D, opt_decomp, para_132D)

%% get response profiles for components
R = Decomp_102D.Rs{1};
plot_on = 1;
% [I_inorder, R_inorder, tags_inorder, snames_inorder] = getResponseProfile(R,plot_on);
[I_inorder, R_inorder, tags_inorder, snames_inorder] = getResponseProfile_ModelMatch(R, plot_on); % need modification

% [I_inorder, R_inorder, tags_inorder, snames_inorder] = getResponseProfile_LZVoc_mix(R,plot_on);
% [I_inorder, R_inorder, tags_inorder, snames_inorder] = getResponseProfile_LZVoc_major(R, plot_on);
% [I_inorder, R_inorder, tags_inorder, snames_inorder] = getResponseProfile_NatVoc(R, plot_on);
% getResponseProfile_ModelMatch(R,plot_on);


%% ========== Misc ============= 
%% View X in 2d pictures
opt.ampLimit    = 0.15.*[-1 1];

figure, h = imagesc(reshape(X(1,:), para.height, para.width), [opt.ampLimit(1), opt.ampLimit(2)]); 
axis image
colorbar
% for i = 1:size(X,1)
for i = 1:180
    set(h,'CData',reshape(X(i,:), para.height, para.width))
    title(num2str(i))
    pause 
end

%% select sound to plot cochleogram & corresponding response pattern
addpath('D:\SynologyDrive\=code=\Sound_analysis\')
load('D:\SynologyDrive\=data=\F_halfcosine_marm_NatJM.mat')
folder_sound = 'D:\SynologyDrive\=sounds=\Vocalization\LZ_AudFilt\MX\';
load('D:\SynologyDrive\=data=\SpecTempParameters_Yueqi')
list = dir(fullfile(folder_sound,'*.wav'));
names_sound = natsortfiles({list.name})';
nSound = 6;

% ind = 1:nSound; % desending order
ind = length(list): -1 : length(list) - nSound + 1; % reversed order

for iComp = 2
    figurex;
%     ('Name',['Component ',num2str(iComp)]), se
t(gcf, 'color','w')
    set(gcf,'Position',[1 41 3440 1323])
    ha = tight_subplot(2, nSound, 0.01, [.1 .01], [0.05 .01]);
    for i = 1:nSound
        % calculate cochleogram from sounds
%         Sd.SoundName = names_sound{I_inorder(ind(i),iComp)};
%         filename = [folder_sound,Sd.SoundName];
%         [Sd.wav,Sd.fs] = audioread(filename);
% %         subplot(1,10,iSound)
%         windur = 0.0025;
%         mode = 'ERB'; % log or linear, or ERB scale
%         plotON = 1;
%         axes(ha(i));
%         %     [F.CochEnv, F.CochEnv_ds, F.CochEnv_dB, F.cf, F.t_ds]  =  getCochleogram(Sd, windur, mode, plotON);
%         [Mat_env_ds, Mat_env, P] =  getCochleogram_halfcosine(Sd, P, plotON); % Mat_env_ds is the compressed cochleagram
%         axis square
%         if i == 1
%         else
%             axis off
%         end
%         drawnow
       
        % get cochleogram directly from F
        axes(ha(i));
        img = F.CochEnv_ds_log(:,:,I_inorder(ind(i), iComp));
        imagesc(flipud(img))
        axis off
        axis square
       
%         xlabel('time(s)'), ylabel('frequency (Hz)')
        title(F.sound_names{I_inorder(ind(i), iComp)},  'Interpreter', 'none')
        
    % plot the response patterns 
%     figure('Name',['Component ',num2str(iComp)]), set(gcf,'color','white')
    opt.ampLimit    = 0.18.*[-1 1];
    opt.trials      = I_inorder(ind(i), iComp);
    img_temp = reshape( X(I_inorder(ind(i), iComp),:), para.height, para.width );
    axes(ha(i+nSound));
    if strcmp(para.side,'RIGHT')
        img_temp = fliplr(img_temp);
    end
    imagesc(img_temp, opt.ampLimit), colormap(jet), axis image, axis off
    colorbar
    drawnow
    end
end

%% visualize X hats (reconstructed X)
i = 4;
figure, 
subplot(1,3,1), imagesc(reshape(X(i, :), para.height, para.width)), axis image
subplot(1,3,2), imagesc(reshape(X_mask(i, :), para.height, para.width)), axis image
subplot(1,3,3), imagesc(reshape(-X_hats{7}(i, :), para.height, para.width)), axis image

% for k = 1:8
%     recon_error_mask(k) = norm((1-M).*(X - X_hats{k}), 'fro');
% end
% figure, plot(recon_error_mask)
%% ============== Generate a mask for component analysis ==============
%======== masks drawn manually ========
% [file_tif,path_tif] = uigetfile('*.tif','Select a surface image','X:\');
% I = imread(fullfile(path_tif,file_tif));
% I = double(imresize(I,para.height./size(I,1)));
I = squeeze(mean(DataMat(:,1,:,:,1),1));
figure, imshow(I,[])
% h = images.roi.Circle(gca,'Center',[floor(para.width/2) floor(para.height/2)],'Radius',floor(para.height/2)); 
% h = images.roi.Rectangle(gca,'Position',[floor(para.width/2)-floor(para.height/2), 1,...
%                                          floor(para.height), floor(para.height)]); 
% h = images.roi.Ellipse(gca,'Center',[floor(para.width/2) floor(para.height/2)],'Semiaxes',[40 20]); 
h = images.roi.Polygon(gca,'Position',[1 1; 1 para.height - 30; 30 para.height; para.width para.height; para.width 1]); 
title('Press Enter after the position is adjusted')
pause
mask_manual = createMask(h);

%======== masks based on response properties ========
% p_th    = 0.0001; % significance threshold
% r_th    = 0.3; % reliability threshold
% c_th    = 0.3; % correlation (across reps) threshold
% plot_on = 1; % plot: 0 = no plots; 1 = plot masks; 2 = plot all figures
% [mask_sig, mask_consis, mask_corr] = PixSelect(para,DataMat,I,p_th,r_th,c_th,plot_on);
mask_final = mask_manual;
% ==== generate final masks with surface image ====
mask_5d = repmat(mask_final,1, 1, para.nRep, para.nStim, para.nFrame);
mask_5d = permute(mask_5d,[3,4,1,2,5]);
DataMat_mask = DataMat.*mask_5d;

%% get response ranking based on applied core- mask
R_core = mean(X,2);
[~, ind_core] = sort(R_core);
% select sound to plot cochleogram & corresponding response pattern
addpath('D:\=code=\Sound_analysis')
folder_sound = 'D:\=code=\McdermottLab\sound_natural\';
list = dir(fullfile(folder_sound,'*.wav'));
names_sound = natsortfiles({list.name})';

figure('color','w')
set(gcf,'Position',[1 41 3440 1323])
ha = tight_subplot(1,10,0.01,[.1 .01],[0.05 .01]);
for i = 1:10
    Sd.SoundName = names_sound{ind_ratio(i)};
    filename = [folder_sound,Sd.SoundName];
    [Sd.wav,Sd.fs] = audioread(filename);
%         subplot(1,10,iSound)
    windur = 0.0025;
    mode = 'ERB'; % log or linear, or ERB scale
    plotON = 1;
    axes(ha(i));
    %     [F.CochEnv, F.CochEnv_ds, F.CochEnv_dB, F.cf, F.t_ds]  =  getCochleogram(Sd, windur, mode, plotON);
    [Mat_env, Mat_env_ds, MatdB, cf, t_ds]  =  getCochleogram(Sd, windur, mode, plotON); % Mat_env_ds is the compressed cochleagram
    axis square
    if i == 1
    else
        axis off
    end
    drawnow
end

% plot the response patterns 
figure('color','white')
opt.ampLimit    = 0.4.*[0 1];
opt.trials      = ind_ratio(1:10);
opt.p           = [1 10];
%     opt.tWindow     = [para.preStim, para.preStim + 12]; % start and end of integration window for calculating response amplitude
[~, ~] = ViewData(DataMat, para, opt); % X may contain NaNs if there are masked pixels



%% plot temporal trace
ampLimit = [-2, 35];
% Stim = [1, 15, 30, 50];
Stim = 129; 

[p,~] = numSubplots(length(Stim));
figurex([1888         302         766         304]);
for i = 1:length(Stim)
    iStim = Stim(i);
%     iStim = 50;
    subplot(p(1), p(2), i)
    Data_temp       = 100.*reshape( DataMat_norm(iStim,:,:,:), para.height*para.width, para.nFrame );
    t = 1/para.fr: 1/para.fr: (para.preStim + para.durStim + para.postStim);
%     ampLimit = [-5, max(Data_temp(:) + 5)];

    hold on, 
    area([para.preStim, para.preStim + para.durStim + 4],...
        [ampLimit(1), ampLimit(1)], ...
        'LineStyle',    ':',...
        'FaceColor',    0.8.*[1 1 1],...
        'FaceAlpha',    0.5);
    area([para.preStim, para.preStim + para.durStim],...
        [ampLimit(2), ampLimit(2)], ...
        'LineStyle',    ':',...
        'FaceColor',    0.8.*[1 1 1],...
        'FaceAlpha',    0.5);
    plot(t, Data_temp')
%     title(['Stim #', num2str(Stim(i))])
    set(gca,'Ylim',ampLimit,'Xlim',[0 max(t)],'fontsize',24)
    set(gca, 'Ytick', [0:10:ampLimit(2)])
    set(gca,'YtickLabels', strcat( string(get(gca,'Ytick')),'%'))
end


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
%% select trials to compare Calcium and intrinsic imaging (temp)
% the entire set of natural sounds
folder_sound = 'D:\=code=\McdermottLab\sound_natural\';
list = dir(fullfile(folder_sound,'*.wav'));
names_sound = natsortfiles({list.name})';
% the subset of natural sounds for intrinsic imaging
folder_sound = 'X:\=Sounds=\Natural_XINTRINSIC2';
list = dir(fullfile(folder_sound,'*.wav'));
names_sound2 = natsortfiles({list.name})';
% pick out trials for the selected natural sounds
opt.trials = [];
for i = 1:length(names_sound2)-1
    opt.trials(i) = find(strcmp(names_sound2(i+1),names_sound));
end

%% frame-to-frame correlation, pick out trials with movement
tic, DataMat = permute(DataMat,[3,4,5,2,1]); time.permute = toc % DataMat = [rep, trial, height, width, frams]
for i = 1:6
tic, DataMat(:,:,:,:,i) = DataMat(:,:,:,para.order(i,:),i); time.reorder = toc % re-arrange according to the experiment order
end
Y = reshape(DataMat, para.height, para.width, para.nFrame*para.nStim*para.nRep); 
Y = Y - min(Y);
nnY = quantile(Y(:),0.005);
mmY = quantile(Y(:),0.995);
[cY,mY,~] = motion_metrics(Y,5); % cY: correlation coefficients between frames

% ==== do registration =======

% ============================
[cM1,mM1,~] = motion_metrics(Yreg,5);

t = 1:para.nFrame*para.nStim*para.nRep;
figure,
ax1 = subplot(2,2,1:2); plot(t, cY), xlim([0, max(t)]), ylim([min(cY), 1]), title('Original')
ax2 = subplot(2,2,3:4); plot(t, cM1), xlim([0, max(t)]), ylim([min(cY), 1]), title('Motion corrected')
linkaxes([ax1,ax2],'xy')
% subplot(2,2,3), histogram(cY, 100)
% subplot(2,2,4), histogram(cY, 100), ylim([0 500])

DataMat_reg = reshape(Yreg, para.height, para.width, para.nFrame, para.nStim, para.nRep);
for i = 1:6
[~,ind] = sort(para.order(i,:));
tic, DataMat_reg(:,:,:,:,i) = DataMat_reg(:,:,:,ind,i); tReorder = toc % re-arrange according to the experiment order
end
tic, DataMat_reg = permute(DataMat_reg,[5, 4, 1, 2, 3]); tPermute = toc % DataMat_reg = [height, width, frames, trial, rep]
