% this script was used for separating core/non-core, and examine response
% profiles for natural sounds

%% get significantly activated area outline
% load data in script_TrialBased first
[DataMat, para] = getDataMat;
opt             = struct;
opt.trials      = 18+[1:9]; % the 60dB trials in tonotopy session 7/26
[X1, DataMat_norm1, DataMat_norm_sep] = getX(DataMat, para, opt);
% opt.p = [1 9];
% opt.ampLimit = [-0.15 0.15];
% figure
% [X1, DataMat_norm1] = ViewData(DataMat, para, opt);
%% ========= variance map (selectivity) ============
map_sel         = selectivity(DataMat_norm1,[30 99]);
options.title   = 'Activation map for tones';
options.method  = 'var';
% options.X_exp = X_exp;
X_tone          = getRespMap(DataMat_norm1, para, options);
%% get response profile for natural sounds
% load data in script_TrialBased first
[DataMat, para] = getDataMat;
opt             = struct;
[X2, DataMat_norm2, ~] = getX(DataMat, para, opt);
%% ========= variance map (across stim, selectivity) ============
map_sel1        = selectivity(DataMat_norm2,[30 98]);
options.title   = 'Activation map for natural sounds';
options.method  = 'var';
X_exp           = getRespMap(DataMat_norm2, para, options);
%% define core using significance measurements
mask_sig = getSigMap(DataMat_norm1, X_tone, para);
xstart = 10; ystart = 11;
mask1 = mask_sig(xstart:end, ystart:end);
% define a noise map to be excluded from non-core
mask_noise = getNoiseMap(DataMat_norm1, X_exp, para);
mask_final = (~mask_noise & ~mask1);

%% load tonotopic and tested trial surface images, get registration parameters
filename_exp    = 'X:\2018.03 T2 (Marmoset 80Z, Xintrinsic, Green & Fluo)\M80Z-180724-13\180724T162507_Blue_Koehler_Fluo_GFP.tif';
filename_tone   = 'X:\2018.03 T2 (Marmoset 80Z, Xintrinsic, Green & Fluo)\M80Z-180726-14\180726T143800_Blue_Koehler_Fluo_GFP.tif';

regopt = struct; 
regopt.mode = 'auto'; regopt.manual_method = 'polynomial'; regopt.auto_method = 'multimodal';
% [img_tone_reg, img_exp, tform] = RegisterSurface(filename_tone, filename_exp, para, regopt);
[img_tone_reg, img_exp, tform] = RegisterSurface(para, regopt);

% img_tone        = imread(filename_tone);
% img_tone        = double(imresize(img_tone,para.height./size(img_tone,1)));
% img_exp         = imread(filename_exp);
% img_exp         = double(imresize(img_exp,para.height./size(img_exp,1)));
% [optimizer,metric]  = imregconfig('multimodal');
% tform               = imregtform(img_tone,img_exp,'rigid', optimizer, metric);
% img_tone_reg        = imwarp(img_tone, tform,'OutputView',imref2d(size(img_exp)));
% figurex;
% subplot(1,2,1), imshowpair(img_exp, img_tone,'Scaling','joint'), title('before registration')
% subplot(1,2,2), imshowpair(img_exp, img_tone_reg,'Scaling','joint'), title('after registration')
% X_tone_reg          = imwarp(X_tone, tform,'OutputView',imref2d(size(img_exp)),'FillValue',0);


%% cut off registration edges
X_diff = log10(X_exp)-log10(X_tone_reg); 
xstart = 10; ystart = 11;
X_diff = X_diff(xstart:end, ystart:end);
X_tone_reg = X_tone_reg(xstart:end, ystart:end);
X_exp = X_exp(xstart:end, ystart:end);
newsize = size(X_diff);

figurex([938         918        2268         420]);
subplot(1,3,1), imagesc(log10(X_tone_reg)), axis image, colormap(jet)
h = colorbar; title(h,'          10^')
title('Activation pattern for tones')
subplot(1,3,2), imagesc(log10(X_exp)), axis image, colormap(jet) 
h = colorbar; title(h,'          10^')
title('Activation pattern for natural sounds')
subplot(1,3,3), imagesc(X_diff), axis image, colormap(jet) 
h = colorbar; title(h,'          10^')
title('log(natural)-log(tone)')


%% generate a mask for "tone-responsive region"
figurex;
mask1 = X_tone_reg > prctile(X_tone_reg(:),89);
mask_outline1 = boundarymask(mask1,4);
imagesc(log10(X_tone_reg.*~mask_outline1)), axis image, colormap(jet)
h = colorbar; title(h,'          10^')
t = text(1,3,['percentile = ', num2str(89)]); t.Color = 'white'; t.FontSize = 14;
mask1 = ~mask1;
%% generate a mask for "non-core"
% figure
% for i = 90:100
i = 90;
mask2 = X_diff > prctile(X_diff(:),i);
mask_outline2 = boundarymask(mask2,4);
figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w');
imagesc(X_diff.*~mask_outline2), axis image, colormap(jet)
h = colorbar; title(h,'          10^')
t = text(1,3,['percentile = ', num2str(i)]); t.Color = 'white'; t.FontSize = 14;
% pause
% end
%% plot final masks with activation pattern
% mask_final = mask1 & mask2 & mask3;
% mask_final = mask1 & mask2;
mask_outline = boundarymask(mask_final,4);

mask_outline1 = double(~boundarymask(mask1,4));
mask_outline1(mask_outline1==0) = 1e7;
Max = max(log10(X_exp(:))); Min = min(log10(X_exp(:))); 
figurex;
imagesc(log10(X_exp.*~mask_outline.*mask_outline1),[Min Max]),axis image, colormap(jet), colorbar
% imagesc(log10(X_exp.*mask_outline1),[Min Max]),axis image, colormap(jet)
h = colorbar; title(h,'          10^')
title('Red: core; Blue: non-core')

MeanDiff_noncore = 10^ (sum(X_diff.*mask_final,'all')/sum(mask_final,'all'));
% MeanDiff_noncore = 10^ (sum(X_diff.*mask1,'all')/sum(mask1,'all'));
MeanDiff_core = 10^ (sum(X_diff.*~mask1,'all')/sum(~mask1,'all'));
t = text(1,5,[' MeanDiff, core = ', num2str(MeanDiff_core),...
    '\newline MeanDiff, noncore = ', num2str(MeanDiff_noncore)]); 
t.Color = 'white'; t.FontSize = 12;


%% generate response profile for masked & unmasked areas
X2_temp = reshape(X2,size(X2,1),para.height, para.width);
X2_cut = X2_temp(:,xstart:end,ystart:end);
X2_cut = reshape(X2_cut,size(X2,1),newsize(1)*newsize(2));
X1_temp = reshape(X1,size(X1,1),para.height, para.width);
X1_cut = X1_temp(:,xstart:end,ystart:end);
X1_cut = reshape(X1_cut,size(X1,1),newsize(1)*newsize(2));

mask_core = reshape(~mask1, [size(mask1,1)*size(mask1,2), 1]);
% mask_belt = reshape(mask_final, [size(mask_final,1)*size(mask_final,2), 1]); % non-core as natural sounds>pure tone
mask_belt = reshape(mask1, [size(mask1,1)*size(mask1,2), 1]); % non-core as everything out of coure
R = X2_cut*mask_belt./length(find(mask_belt)) ./ (X2_cut*mask_core./length(find(mask_core)));
% X = X2_cut;
% R = [];
% R(:,3) = sum(X,2);
% % R(:,1) = X*mask_core./length(find(mask_core)); % averaged response, core
% R(:,1) = 100.*X*mask_belt./length(find(mask_belt)) ./ (X*mask_core./length(find(mask_core)));
% R(:,2) = X*mask_belt./length(find(mask_belt));
% 
% figurex([1440 804 680 534]);
% hold on, 
% scatter(R(:,1), R(:,2))
% [p,rsq,~] = RSquared(R(:,1),R(:,2))
% % xlabel('Mean response amplitude in core')
% xlabel('Relative response in non-core (%)')
% ylabel('Mean response amplitude in non-core')
% axis square
% % ylim(get(gca,'xlim'))
% % 
% R = [];
% X = X1_cut;
% R(:,3) = sum(X,2);
% % R(:,1) = X*mask_core./length(find(mask_core)); % averaged response, core
% R(:,1) = 100.*X*mask_belt./length(find(mask_belt)) ./ (X*mask_core./length(find(mask_core)));
% R(:,2) = X*mask_belt./length(find(mask_belt));
% scatter(R(:,1), R(:,2))
% [p,rsq,~] = RSquared(R(:,1),R(:,2))
% legend({'Natural sound session', 'Tone session'})

% R = (X2*~mask_linear./length(find(~mask_linear)))...
%     ./(X2*mask_linear./length(find(mask_linear)));

%%
[~, ind_ratio] = sort(R(:,5));
% select sound to plot cochleogram & corresponding response pattern
addpath('D:\=code=\Sound_analysis')
folder_sound = 'D:\=code=\McdermottLab\sound_natural\';
list = dir(fullfile(folder_sound,'*.wav'));
names_sound = natsortfiles({list.name})';

figure('color','w','Position',[1 41 3440 800])
ha = tight_subplot(2,10,0.01,[.05 .01],[0.05 .01]);
for i = 1:20
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
figurex;
opt = struct;
opt.ampLimit    = 0.4.*[-1 1];
opt.trials      = ind_ratio(1:20);
opt.p           = [2 10];
opt.saveON = 0;
%     opt.tWindow     = [para.preStim, para.preStim + 12]; % start and end of integration window for calculating response amplitude
[~, ~] = ViewData(DataMat, para, opt); % X may contain NaNs if there are masked pixels

%% project tone-responsive areas to components
% for i = 1:9
% X_temp = reshape(X1(i,:), para.height, para.width);
% X_temp = imwarp(X_temp, tform,'OutputView',imref2d(size(img_exp)),'FillValue',0);
% figure, imagesc(X_temp);
% X1_reg(i,:) = reshape(X_temp,1,para.height*para.width);
% end

ind_comp = [3 2 4 6 5 1]; % for 80Z
figure,set(gcf,'color','white')
freqs = [0.11, 0.22, 0.44, 0.88, 1.76, 3.52, 7.04, 14.08, 28.16];
% figure, set(gcf,'color','white')
ind = 1:9;
for k = 1:length(ind_comp)

    for i = 1:9
%     subplot(1,9,i)
%     imagesc(reshape(X1(ind(i),:), para.height, para.width), [0 0.2]); axis image, colormap(jet)
        temp = corrcoef(W(ind_comp(k),:),X1_reg(ind(i),:));
        proj(i) = temp(1,2);
    end
    
    subplot(1,6,k)
    plot(freqs, proj, 'linewidth',4,'Marker','x')    
    hold on, plot([freqs(1), freqs(end)], [0 0],'linestyle','--','color','k')    

    set(gca,'xscale','log');
    ylim([-0.6 0.8]);
    xlim([freqs(1) freqs(end)]);
    set(gca,'xtick', freqs(1:2:end))
    set(gca,'ytick',-0.6:0.2:0.6)
    set(gca,'fontsize',24);
    axis square
    xlabel('Frequency (kHz)','fontsize',24),
    %     title(['Correlation coefficient, component #',num2str(ind(iComp))],'fontsize',14)
    xtickangle(45)
end


%% regression with sound features
% load('D:\=code=\Sound_analysis\F_test.mat') 
K = 1;
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
    xlabel('Frequency (Hz)','fontsize',24),
%     title(['Correlation coefficient, component #',num2str(ind(iComp))],'fontsize',14)
    xtickangle(45)

    subplot(2,K,i + K), 
%     spectemp_r = reshape(result_r(F.nFreq+1:F.nFreq+F.nSpectemp,ind(i)), ...
%                         size(F.spectemp_mod,1), size(F.spectemp_mod,2));
    spectemp_r = reshape(result_r(F.nFreq+F.nTemp+F.nSpec+1:F.nFreq++F.nTemp+F.nSpec+F.nSpectemp,ind(i)), ...
                        size(F.spectemp_mod,1), size(F.spectemp_mod,2));
    cmax = max(max(  abs( result_r(F.nFreq+1 : F.nFreq+F.nSpectemp,1:K) )  )); 
    imagesc(flipud(spectemp_r),[-cmax, cmax]), colormap('jet')
    colorbar
    set(gca,'ytick',1:2:length(F.spec_mod_rates))
    set(gca,'yticklabels',arrayfun(@num2str,fliplr(F.spec_mod_rates(1:2:end)),'UniformOutput',false))
    set(gca,'xtick',1:2:length(F.temp_mod_rates)-1)
    set(gca,'xticklabels',arrayfun(@num2str,F.temp_mod_rates(2:2:end),'UniformOutput',false))
    axis square
    set(gca,'fontsize',24);
end
%%
figurex([1440         484         560         854])
subplot(2,1,1)
plot(F.temp_mod_rates(2:end), result_r(F.nFreq+1:F.nFreq+F.nTemp,ind(i)), 'linewidth',4,'Marker','x')
    hold on, plot([F.temp_mod_rates(2), F.temp_mod_rates(end)], [0 0],'linestyle','--','color','k')    
    set(gca,'xscale','log');
    ymax = max(max(  abs( result_r(F.nFreq+1:F.nFreq+F.nTemp,1:K) )  )); ylim([-0.77,0.77]);
    set(gca,'xtick',F.temp_mod_rates(2:2:end))
    set(gca,'ytick',-0.6:0.2:0.6)
    set(gca,'xticklabels',arrayfun(@num2str,F.temp_mod_rates(2:2:end),'UniformOutput',false))
    set(gca,'fontsize',24);
    axis square
    xlabel('Temporal modulation rates (Hz)','fontsize',24),
%     title(['Correlation coefficient, component #',num2str(ind(iComp))],'fontsize',14)
subplot(2,1,2)
plot(1:F.nSpec, result_r(F.nFreq+F.nTemp+1:F.nFreq+F.nTemp+F.nSpec,ind(i)), 'linewidth',4,'Marker','x')
    hold on, plot([F.spec_mod_rates(1), F.spec_mod_rates(end)], [0 0],'linestyle','--','color','k')    
%     set(gca,'xscale','log');
    ylim([-0.77,0.77]);
    xlim([1, F.nSpec])
    set(gca,'xtick',1:2:F.nSpec)
    set(gca,'ytick',-0.6:0.2:0.6)
    set(gca,'xticklabels',arrayfun(@num2str,F.spec_mod_rates(1:2:end),'UniformOutput',false))
    set(gca,'fontsize',24);
    axis square
    xlabel('Spectral modulation rates (Cyc/Oct)','fontsize',24),
%     title(['Correlation coefficient, component #',num2str(ind(iComp))],'fontsize',14)
%%
function map = selectivity(DataMat_norm,prctilerange)
    temp = squeeze(mean(DataMat_norm,4));
    map_mean = squeeze(mean(temp,1));
    map_var = squeeze(var(temp,[],1));
    map = map_var./map_mean; map(map<0) = NaN; map = log10(map);
    Max = prctile(map(:), prctilerange(2));
    Min = prctile(map(:), prctilerange(1));
    figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w');
    imagesc(map, [Min, Max]), axis image; 
    h = colorbar; title(h,'          10^')
    title('Var. of response pattern')

end

function X = getRespMap(DataMat_norm, para, options)
if strcmp(options.method, 'mean')
    % ========= use mean amplitude as a measurement ====
    X = reshape(mean(X1,1), para.height, para.width);
    figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w');
    imagesc(X), axis image, colormap(jet)
    title(options.title)
elseif strcmp(options.method, 'var')
    % ==== use variance as a measurement ====
    temp = permute(DataMat_norm,[2 3 4 1]);
    temp = reshape(temp,[para.height, para.width, size(DataMat_norm,1)*para.nFrame]);
    X = var(temp,[],3);
    figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w');
    if isfield(options,'X_exp')
        Max = max(log10(options.X_exp(:)));Min = min(log10(options.X_exp(:)));
    else
        Max = max(log10(X(:)));Min = min(log10(X(:)));
    end
    imagesc(log10(X), [Min Max]), colormap(jet), axis image; 
    h = colorbar; title(h,'          10^')
    title(options.title)
else
    disp('set 3rd variable as mean or var')
end
end

function mask = getSigMap(DataMat_norm, X_tone, para)
% get significance map
% =========== mask based on spatial distribution ========
figure, imagesc(log10(X_tone)), axis image, colormap(jet)
h = images.roi.Rectangle(gca,'Position',[0.5000 0.5000 32.4885 75]); 
title('Press Enter after the position is adjusted')
pause
mask_manual = createMask(h);
temp1 = mask_manual.*X_tone;
mean1 = mean(temp1,'all');
std1 = std(temp1,[],'all');
mask = X_tone > 50*mean1;
mask_outline = boundarymask(mask,4);
figurex, imagesc(log10(X_tone.*~mask_outline)), axis image, colormap(jet), colorbar
% =========== mask based on temporal trace ========
% temp = permute(DataMat_norm,[2 3 4 1]);
% % temp1 = [reshape(temp(:,:,1:para.fr,:), para.height*para.width, length(opt.trials)*para.fr),...
% %     reshape(temp(:,:,para.fr*3+1:para.fr*5,:), para.height*para.width, length(opt.trials)*2*para.fr)];
% temp1 = reshape(temp(:,:,1:para.fr,:), para.height*para.width, size(temp,4)*para.fr);
% temp2 = reshape(temp(:,:,para.fr+1:para.fr*3,:), para.height*para.width, size(temp,4)*2*para.fr);
% mean2 = mean(temp2,2);
% mean1 = mean(temp1,2); 
% std1 = std(temp1,[],2);
% sig_map = mean2 > mean1+5*std1;
% mask = reshape(sig_map, para.height, para.width);
% mask_outline = boundarymask(mask,4);
% figurex, imagesc(log10(X_tone.*~mask_outline)), axis image, colormap(jet), colorbar
end

function mask = getNoiseMap(DataMat_norm, X, para)
% get noise map to be excluded from non-core
% =========== mask based on spatial distribution ========
figure, imagesc(log10(X)), axis image, colormap(jet)
h = images.roi.Rectangle(gca,'Position',[0.5000 0.5000 24.7120 75]); 
title('Press Enter after the position is adjusted')
pause
mask_manual = createMask(h);
temp1 = mask_manual.*X;
mean1 = mean(temp1,'all');
std1 = std(temp1,[],'all');
mask = X < 50*mean1;
mask_outline = boundarymask(mask,4);
figurex, imagesc(log10(X.*~mask_outline)), axis image, colormap(jet), colorbar
end