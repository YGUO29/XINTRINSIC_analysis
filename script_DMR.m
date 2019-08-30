clear all

% load recording session
[file_mat,path_mat] = uigetfile('*.mat','Select a mat file to analyze','\\FANTASIA-DS3617\Test_Imaging\2018.12.T1 (Marmoset 96B, Xintrinsic)');
load(fullfile(path_mat,file_mat));clc
[~,nametemp,~,] = fileparts(file_mat);
Sys = load(fullfile(path_mat,[nametemp(1:end-3),'.mat']));


% I = double(imresize(I,1/16));
para.nRep =     size(P.ProcDataMat,1);
para.nStim =    size(P.ProcDataMat,2);
para.height =   size(P.ProcDataMat,3);
para.width =    size(P.ProcDataMat,4);
para.nFrame =   size(P.ProcDataMat,5);


para.fr =       P.ProcFrameRate; % frame rate
para.preStim =  Sys.S.TrlDurPreStim;
para.durStim =  Sys.S.TrlDurStim;
para.postStim = Sys.S.TrlDurPostStim;
para.filename = file_mat;
para.pathname = path_mat;
DataMat =       P.ProcDataMat;

% load stimulus information
[file_mat,path_mat] = uigetfile('*.mat','Select a mat file to analyze','D:\=sounds=\DMR');
load(fullfile(path_mat,file_mat));

if contains(lower(Sys.S.SesSoundFile), 'down')
    para.direction = 'down';
else
    para.direction = 'up';
end

para.durStim = S.dur; 
para.preStim = 0; 
para.postStim = 0;
clear Sys P
size(DataMat)



%% View data
opt.ampLimit    = [-0.15, 0.15];
opt.saveON      = 0; 
opt.soundON     = 0;
ViewData_raw(DataMat, para, opt);
        
%% cycle-based tonotopy analysis

% ===== select frequency component to analysis =====
% period          = S.tm_period; 
% rep             = S.tm_cycles;

period          = para.preStim + para.durStim + para.postStim;
rep             = para.nRep;
% ==================================================


% change the polarity here!!
data_temp       = - squeeze(DataMat); 
if para.nRep == 1
    data_temp   = permute(data_temp,[3 1 2]);
else
    data_temp   = permute(data_temp,[4 1 2 3]);
end
nPix            = para.height*para.width;
data_temp       = reshape(data_temp,[para.nRep*para.nFrame, para.height*para.width]);

% FFT
L               = floor(para.fr * para.nRep * (para.preStim + para.durStim + para.postStim)); %length of signal
T               = 1/para.fr; % sampling period
t               = (0:L-1)*T;

% ========== plot average temporal trace ========== 
data_temp_mean  = mean(data_temp,2);
Max             = max(data_temp_mean) + 3*std(data_temp_mean);
Min             = 2*mean(data_temp_mean) - Max; 
% data_temp_mean  = data_temp(: , sub2ind([para.height, para.width],40,61));
ind             = para.fr.*[period: period: floor(period.*rep)];
ind             = [1,ind];
 
fig = figure; size_scr = get(0,'ScreenSize'); set(gcf,'position',[1 1 size_scr(3:4)])
subplot(2,3,1:3)
hold on
for i = 1:length(ind)-1
    if mod(i,2)
        plot( t(ind(i):ind(i+1)), data_temp_mean(ind(i):ind(i+1)), 'k');
    else
        plot( t(ind(i):ind(i+1)), data_temp_mean(ind(i):ind(i+1)), 'r');
    end
end
set(gca,'XTick', period: period: floor(period.*rep))
% set(gca,'XTick', unique(sort([S.fm_period:S.fm_period:para.durStim, S.tm_period:S.tm_period:para.durStim])) )
% plot(1/para.fr:1/para.fr:para.durStim,  data_temp_mean);
grid on
% =================================================

if mod(L,2)
    
    
    data_temp   = [data_temp;zeros(1,nPix)];
    L = L +1;
end
data_fft                    = fft(data_temp, L, 1);
data_fftamp_mat             = abs( data_fft/L );
data_fftamp_mat             = data_fftamp_mat(1:L/2+1,:);
data_fftamp_mat(2:end-1,:,:)= 2*data_fftamp_mat(2:end-1,:);
data_fftamp_mat             = data_fftamp_mat./repmat(data_fftamp_mat(1,:),[floor(L/2)+1,1]); % normalized to mean amplitude
data_fftamp_mean            = mean(data_fftamp_mat,2);
f                           = para.fr*(0:(L/2))/L;

% ========== plot average spectrum ==========
    subplot(2,4,5)
    semilogx(f(2:end),data_fftamp_mean(2:end)), xlim([0 1])
    % label fm and tm periods
    ind_fm = floor(interp1(f,1:length(f),1/S.fm_period));
    ind_tm = floor(interp1(f,1:length(f),1/S.tm_period));
    hold on, scatter(f([ind_fm, ind_tm]), data_fftamp_mean([ind_fm, ind_tm]))
    title('Averaged spectrum')
% ===========================================

freq_comp = floor(interp1(f,1:length(f),1/period));
% data_fftagl_mat             = angle(data_fft')';
data_fftagl_mat             = mod( angle(data_fft')', 2*pi );
% data_fftagl_mat             = angle(data_fft(1:L/2+1,:,:));
% data_fftagl_mat             = data_fftagl_mat + pi.*ones(size(data_fftagl_mat));

% take the corresponding frequency component
data_fftamp_one = squeeze(data_fftamp_mat(freq_comp,:));
data_fftagl_one = squeeze(data_fftagl_mat(freq_comp,:));
% data_fftamp_one = reshape(data_fftamp_one, nPix, 1);
% data_fftagl_one = reshape(data_fftagl_one, nPix, 1);

% ========== plot phase ====================
    subplot(2,4,6) 
    hist(data_fftagl_one./pi,100);
    title('Phase distribution (uncompensated)')
% ===========================================

% Compensate the pseudo-delay
delay   = 2.6;
hue     = mod(data_fftagl_one - (delay/period)*2*pi, 2*pi)./(2*pi);

% Reverse the hue for DOWN cycle
if strcmp(para.direction,'down')
    hue = 1 - hue;
end
hue_lim = 0.8;
hue = hue.*hue_lim;

saturation      = data_fftamp_one;
% ========== plot amplitude map ================
    h1 = subplot(2,4,7);
    map_amp = reshape(saturation, para.height, para.width); 
    imagesc(map_amp);
    axis image, colormap(h1, gray), colorbar
% ===========================================
saturation_lim  = 0.03; 
saturation      = min(saturation./saturation_lim,1);
map_rgb         = hsv2rgb([hue; saturation; saturation]');
map_rgb         = reshape(map_rgb, [para.height, para.width,3]);

% ========== plot phase map ================
    h2 = subplot(2,4,8);
    imagesc(map_rgb); 
    axis image, colormap(h2, hsv), colorbar
% ===========================================
% % ============== Sing fm =============
%      colorbar('Ticks',0.8.*(para.preStim + para.durStim.*[0:8]./8)./(para.preStim + para.durStim + para.postStim),...
%          'TickLabels',{'4','6.8','8','6.8','4','1.28','0.125','1.28','4'})
% ============== SingLog fm =============
% colorbar('Ticks',0.8.*(para.preStim + para.durStim.*[0:8]./8)./(para.preStim + para.durStim + para.postStim),...
%          'TickLabels',{'1','4.35','8','4.35','1','0.23','0.125','0.23','1'})
 % ============== Sine tm =============
colorbar('Ticks',0.8.*(para.preStim + para.durStim.*[0:8]./8)./(para.preStim + para.durStim + para.postStim),...
         'TickLabels',{'0','22.6','32','22.6','0','-22.6','-32','-22.6','0'})

%% plot two amplitude maps together
map_temp = zeros(para.height, para.width, 3);
saturation_lim = 0.03;
period          = S.fm_period; % select frequency component to analysis
freq_comp = floor(interp1(fig,1:length(fig),1/period));
saturation = squeeze(data_fftamp_mat(freq_comp,:));
saturation = min(saturation./saturation_lim,1);
saturation = reshape(saturation, para.height, para.width); 
map_temp(:,:,1) = saturation;

period          = S.tm_period; % select frequency component to analysis
freq_comp = floor(interp1(fig,1:length(fig),1/period));
saturation = squeeze(data_fftamp_mat(freq_comp,:));
saturation = min(saturation./saturation_lim,1);
saturation = reshape(saturation, para.height, para.width); 

map_temp(:,:,2) = saturation;

figure, imagesc(map_temp), 
axis image
     
     