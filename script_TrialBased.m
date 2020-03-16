% =============== load files to analyze ================
clear all
para.nRep = 0;  
% 80Z Calcium: file_mat = '180724T140401_Blue_Koehler_Fluo_GFP_P1.mat';
% 132D Calcium 1: file_mat = '190119T114434_Blue_Koehler_Fluo_GFP_P1.mat';
% 132D Calcium 2: file_mat = '190224T095545_Blue_Koehler_Fluo_GFP_P1.mat';
%% add more repetitions
[file_mat,path_mat] = uigetfile('*.mat','Select a mat file to analyze','X:\', 'Multiselect', 'on');

if ~iscell(file_mat) % only load one file
        load(fullfile(path_mat,file_mat));
        file_parts = strsplit(file_mat,{'_','.'});
        ind = find(strcmp(file_parts,'P1'));
        nametemp = strjoin(file_parts(1:ind-1),'_');
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
            para.filename = file_mat;
            para.pathname = path_mat;
        else
            DataMat(para.nRep+1 : para.nRep+size(P.ProcDataMat,1),:,:,:,:) = squeeze(P.ProcDataMat);
            para.nRep = para.nRep + 1;
        end
        clear S P
else % load multiple files
    for i = 1:length(file_mat)
        load(fullfile(path_mat,file_mat{i}));
        file_parts = strsplit(file_mat{i},{'_','.'});
        ind = find(strcmp(file_parts,'P1'));
        nametemp = strjoin(file_parts(1:ind-1),'_');
        load(fullfile(path_mat,[nametemp,'.mat']));

    if para.nRep + i == 1 % read the initial file
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
    
    clear S P
    para.nRep = size(DataMat,1);
    disp([num2str(para.nRep),' finished'])
    end     
end

%% check alignment across sessions
test_frames = DataMat(:, [1 20 30], :,:,:);
test_frames = permute(test_frames,[3,4,1,2,5]);
Y = reshape(test_frames, [para.height, para.width, size(test_frames,3)*size(test_frames,4)*size(test_frames,5)]); 
figure,
for i = 1:2000
    imshow(Y(:,:,i),[]),title(num2str(i))
    pause(0.1)
end
Y = Y - min(Y(:));
nnY = quantile(Y(:),0.005);
mmY = quantile(Y(:),0.995);
[cY,mY,~] = motion_metrics(Y,5); % cY: correlation coefficients between frames
t = 1:size(test_frames,3)*size(test_frames,4)*size(test_frames,5);
figure,
plot(t, cY), xlim([0, max(t)]), ylim([min(cY), 1]), title('Original')

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
%% View data
% ind_best = [79 163 88 15 4 155 154 97 83 75];
% ind_worst = [26 89 142 158 100 8 32 123 108 99];

% for i = 1:10
    tic
    opt.ampLimit    = 0.01.*[-1 1];
    opt.mode        =  'avgrep'; % avgrep or allrep
    opt.plotMode    =  'combined'; % combined or separate (video saving is only available for 'combined' mode)
%     opt.trials =    ind; 
%     opt.trials      = 1:para.nStim;
    opt.saveON      = 0; 
    opt.soundON     = 0;
    opt.reps        = [];
    opt.tWindow     = [para.preStim, para.preStim + 18]; % start and end of integration window for calculating response amplitude
%     opt.tWindow     = []; % start and end of integration window for calculating response amplitude
    figure, set(gcf, 'color','w')
    [X, DataMat_norm] = ViewData(DataMat, para, opt); % X may contain NaNs if there are masked pixels
    % X = getX(DataMat, para, opt);
    toc
% end
% ======== construct a X without NaN ========
[~, ind_delete]     = find( isnan(X) ); % linear index
ind_save            = setdiff(1:para.width*para.height, ind_delete);
X(:, ind_delete)    = [];

%% perform ICA analysis

addpath('D:\=code=\McdermottLab\toolbox_nonparametric-ICA') % MAKE SURE
% THIS FOLDER IS ON THE TOP OF THE PATH LIST IN MATLAB!!!!
K = 6;
RANDOM_INITS = 10;      
PLOT_FIGURES = 0;
% ========= reverse the sign for X for intrinsic imaging ========
tic,[R, W] = nonparametric_ica(-X, K, RANDOM_INITS, PLOT_FIGURES);toc

comp = cell(1,K);
% I_norm = (I - min(min(I)))./(max(max(I)) - min(min(I)));

cutoff = 0.005;
% [~,ind] = sort(sum(abs(W),2),'descend');
% ind = [3 2 4 6 5 1]; % for 80Z
% ind = [4 2 3 5 1 6]; % for 132D, session 2
ind = 1:K;
figure, set(gcf, 'color','w')
[p,n] = numSubplots(K);
for i = 1:K
    cutoff = mean(W(i,:)) + 7*std(W(i,:)); % variable cutoff values for each components
    comp{i} = zeros(para.height, para.width);
    comp{i}(ind_save) = W(ind(i),:);
% %     % ============= plot components only =============
    subplot(p(1),p(2),i),imagesc(comp{i},cutoff.*[-1 1]),axis image, colormap(jet)
    colorbar,
    axis off
    % ============= plot components with image =============
%     mask = 0.7.*comp{i}./max(abs(comp{i})); 
%     mask = comp{i};
%     mask(mask > cutoff) = cutoff; mask(mask < - cutoff) = - cutoff;
%     mask = mask.*8; 
%     img = repmat(I_norm,1,1,3); % three layers, representing R,G,B 
%     img(:,:,1) = img(:,:,1) + mask;
%     subplot(p(1),p(2),i),
%     imagesc(img),axis image
%     axis off
end

%% get response profiles for components
plot_on = 1;
[I_inorder, R_inorder, tags_inorder, snames_inorder] = getResponseProfile(R,plot_on);
%%
addpath('D:\=code=\Sound_analysis')
folder_sound = 'D:\=code=\McdermottLab\sound_natural\';
list = dir(fullfile(folder_sound,'*.wav'));
names_sound = natsortfiles({list.name})';
ind = 165-9:165;
for iComp = 1
    figure('Name',['Component ',num2str(iComp)]), set(gcf, 'color','w')
    set(gcf,'Position',[1 41 3440 1323])
    ha = tight_subplot(1,10,0.01,[.1 .01],[0.05 .01]);
    for i = 1:10
        Sd.SoundName = names_sound{I_inorder(ind(i),iComp)};
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
    figure('Name',['Component ',num2str(iComp)]), set(gcf,'color','white')
    opt.ampLimit    = 0.02.*[-1 1];
    opt.mode        =  'avgrep'; % avgrep or allrep
    opt.plotMode    =  'combined'; % combined or separate (video saving is only available for 'combined' mode)
    % opt.trials =    44+[1:12,23:28]; 
%     opt.trials      = I_inorder(1:10, iComp);
    opt.trials      = I_inorder(ind, iComp);
    opt.saveON      = 0; 
    opt.soundON     = 0;
    opt.reps        = [];
    opt.tWindow     = [para.preStim, para.preStim + 12]; % start and end of integration window for calculating response amplitude
    [~, ~] = ViewData(DataMat_mask, para, opt); % X may contain NaNs if there are masked pixels
end

   

%% regression with sound features

load('D:\=code=\Sound_analysis\F_yg_marm.mat') 
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
f = figure; 
set(gcf,'position',[1,1,0.8.*scr_size]); 
set(gcf, 'color','w')
% ind = [1 2 5 6 4 3];
% ind = [1 2 3 4 5 6]; % plot order
% ind = [3 2 4 6 5 1];
ind = 1:K;
for i = 1:K
    iComp = ind(i);
% for iComp = 1:K
    % correlation with frequency power
%     figure(f1)    
    subplot(2,K,i),
    
    plot(F.FreqBounds(1:end-1), result_r(1:F.nFreq,ind(i)), 'linewidth',4,'Marker','x')
    hold on, plot([F.FreqBounds(1:2), F.FreqBounds(end-1)], [0 0 0],'linestyle','--','color','k')    
    set(gca,'xscale','log');
    ymax = max(max(  abs( result_r(1:F.nFreq,1:K) )  )); ylim([-0.77,0.77]);
    set(gca,'xtick',F.FreqBounds(1:2:end))
    set(gca,'ytick',-0.6:0.2:0.6)
    set(gca,'xticklabels',arrayfun(@num2str,F.FreqBounds(1:2:end-1),'UniformOutput',false))
    set(gca,'fontsize',24);
    axis square
%     xlabel('Frequency','fontsize',14),
%     title(['Correlation coefficient, component #',num2str(ind(iComp))],'fontsize',14)
    xtickangle(45)
%     figure(f2)

    subplot(2,K,i + K), 
    
    spectemp_r = reshape(result_r(F.nFreq+1:F.nFreq+F.nSpectemp,ind(i)), ...
                        size(F.spectemp_mod,1), size(F.spectemp_mod,2));
    cmax = max(max(  abs( result_r(F.nFreq+1 : F.nFreq+F.nSpectemp,1:K) )  )); 
    imagesc(flipud(spectemp_r),[-cmax, cmax]), colormap('jet')
%     colorbar
    set(gca,'ytick',1:2:length(F.spec_mod_rates))
    set(gca,'yticklabels',arrayfun(@num2str,fliplr(F.spec_mod_rates(1:2:end)),'UniformOutput',false))
    set(gca,'xtick',1:2:length(F.temp_mod_rates)-1)
    set(gca,'xticklabels',arrayfun(@num2str,F.temp_mod_rates(2:2:end),'UniformOutput',false))
    axis square
    set(gca,'fontsize',24);

%     xlabel('Temporal modulation rate (cycles/s)','fontsize',14),
%     ylabel('Spectral modulation rate (cycles/octave)','fontsize',14),
%     title(['Correlation coefficient, component #',num2str(ind(iComp))],'fontsize',14)

%   correlation with the full spectrotemporal modulation power map
%     subplot(2,K,iComp + K), 
%     spectemp_r = reshape(result_r(F.nFreq+1:F.nFreq+F.nSpectemp_full,ind(iComp)), ...
%                         size(F.spectemp_mod_full,1), size(F.spectemp_mod_full,2));
%     cmax = max(max(  abs( result_r(F.nFreq+1 : F.nFreq+F.nSpectemp_full,1:K) )  )); 
%     imagesc(flipud(spectemp_r),[-cmax, cmax]),colorbar, colormap('jet')
%     set(gca,'yticklabels',arrayfun(@num2str,fliplr(F.spec_mod_rates),'UniformOutput',false), 'fontsize',12)
%     set(gca,'xticklabels',arrayfun(@num2str,F.temp_mod_rates,'UniformOutput',false), 'fontsize',12)
%     xlabel('Temporal modulation rate (cycles/s)','fontsize',14),
%     ylabel('Spectral modulation rate (cycles/octave)','fontsize',14),
%     title(['Correlation coefficient, component #',num2str(ind(iComp))],'fontsize',14)

%     saveas(f,['80Z_session1_component_',num2str(iComp),'.png'])

end
% saveas(f1,'132D_session1_FreqPower.png')
% saveas(f2,'132D_session1_SpecTemp.png')
% saveas(f,'80Z_session1_reg.png')

%% Plot responsive area (variance map)
DataMean = squeeze(mean(DataMat_norm,1));
DataVar = var(DataMean,[],3);
% figure,imagesc(DataVar,[0 0.7*max(DataVar(:))]), axis image
figure,imagesc(log10(DataVar)), axis image
axis off
h = colorbar;
set(h,'ticklabels',num2cell(10.^(get(h,'ticks'))))


%% plot background pictures
figure, set(gcf, 'color','w')
for i = 1:K
    imshow(I_norm,[]),axis image
    axis off
end

% plot a colorbar only
fig1=figure; set(gcf, 'color','w')
axis off
colormap(hsv(100));
caxis([-0.34 0.34]);
h = colorbar('Ticks',[-0.34 0 0.34],...
    'location','Southoutside',...
    'position',[0.2 0.2 0.5 0.1],...
    'FontSize',26);
%     'Ticklabels',{'0%', '45%'},...

%% plot temporal trace
ampLimit = [-5 5];
Stim = 1;

[p,~] = numSubplots(length(Stim));
figure,
for i = 1:length(Stim)
%     iStim = Stim(i);
    iStim = 12;
    subplot(p(1), p(2), i)
    Data_temp       = 100.*reshape( DataMat_norm(iStim,:,:,:), para.height*para.width, para.nFrame );
    t = 1/para.fr: 1/para.fr: (para.preStim + para.durStim + para.postStim);
%     ampLimit = [-5, max(Data_temp(:) + 5)];

    hold on, 
    area([para.preStim, para.preStim + para.durStim],...
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
    set(gca,'Ylim',ampLimit,'Xlim',[0 max(t)],'fontsize',20)
    set(gca,'YtickLabels', strcat( string(get(gca,'Ytick')),'%'))
    set(gcf,'position',[200 200 700 180])
end


%%  Variance across reps 1
opt.plotON      = 1;
opt.ampLimit    = 75; % the percentile of the data distribution as upper limit
opt.tWindow     = []; % select time window for analysis (default is the stimulus-on period)
Var             = AnalysisVar(para,DataMat,opt);

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
