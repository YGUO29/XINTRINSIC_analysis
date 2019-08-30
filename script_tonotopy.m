%% Tonotopy map
% -------------------
% Xindong's code below 
% -------------------
% R.PhPw_CycleAmp =      abs(        R.PhPwQt_FFTRaw(:,:,R.N_Qfc) )*2/R.N_Ft;
% R.PhPw_CycleAgl =      mod(angle(  R.PhPwQt_FFTRaw(:,:,R.N_Qfc) ), 2*pi);                      
% R.PtOne_CycleAmp =     reshape(R.PhPw_CycleAmp,   R.N_Pt, 1);
% R.PtOne_CycleAgl =     reshape(R.PhPw_CycleAgl, R.N_Pt, 1);  
% 
% % Compensate the pseudo-delay
% R.PtOne_Hue =          mod(    R.PtOne_CycleAgl - ...
%                                     T.PsuedoDelay/S.SesSoundDurTotal*2*pi, 2*pi)/(2*pi);
% % Reverse the hue for DOWN cycle
% if contains(lower(S.SesSoundFile), 'down')
%     R.PtOne_Hue =      1 - R.PtOne_Hue;
% end
% 
% % Phase cut    
%     R.PtOne_Hue =      (R.PtOne_Hue-S.TrlDurPreStim/S.TrlDurTotal) /...
%                                 (S.TrlDurStim/S.TrlDurTotal);
%                             % match the Stimulus ONSET / OFFSET to 0-1
%     R.PtOne_Hue =      min(max(R.PtOne_Hue, 0), 1);
%                             % limit the hue range as 0-HueLim
% % CONTINUOUS or DIScontinous hue    
% if contains(lower(S.SesSoundFile),  'sinusoidcycle')
%     R.HueLim =         1.0;
% else
%     R.HueLim =         0.8;
% end 
%     R.PtOne_Hue =      R.PtOne_Hue * R.HueLim;
% 
% % Construct the entire image       
% R.SaturationLim =  0.005;
% R.PtOne_Saturation =   min(R.PtOne_CycleAmp/R.SaturationLim,1);
% R.PtThree_TuneMap =    hsv2rgb([   R.PtOne_Hue,...
%                                         R.PtOne_Saturation,...
%                                         R.PtOne_Saturation]);
% R.PhPwThree_TuneMap =  reshape(R.PtThree_TuneMap, R.N_Ph, R.N_Pw,3);

%%
clear all

MonkeyID = 0; % 1: 80Z; 2: 132D session11; 3: 132D session2; 4: 132D 1st Scrambled; 6: 132D 2nd Scrambled

switch MonkeyID 
    case 1 % 80Z
    file_mat = '180724T133906_Blue_Koehler_Fluo_GFP_P1.mat';
    file_tif = '180724T133546_Blue_Koehler_Fluo_GFP.tif';
    path_mat = '\\FANTASIA-DS3617\Test_Imaging\2018.03 T2 (Marmoset 80Z, Xintrinsic, Green & Fluo)\M80Z-180724-13\';
    load(fullfile(path_mat,file_mat));
    [~,nametemp,~,] = fileparts(file_mat);
    load(fullfile(path_mat,[nametemp(1:end-3),'.mat']));
    I = imread(fullfile(path_mat,file_tif));
    
    case 2 % 132D, session 1, up
    file_mat = '190119T112530_Blue_Koehler_Fluo_GFP_P1.mat';
    file_tif = '190119T132947_Blue_Koehler_Fluo_GFP.tif';
    path_mat = '\\FANTASIA-DS3617\Test_Imaging\2018.11.T1 (Marmoset 132D, Xintrinsic)\M132D-190119-10\';
    load(fullfile(path_mat,file_mat));
    [~,nametemp,~,] = fileparts(file_mat);
    load(fullfile(path_mat,[nametemp(1:end-3),'.mat']));
    I = imread(fullfile(path_mat,file_tif));
    
    case 3 % 132D, session 1, down
    file_mat = '190119T113302_Blue_Koehler_Fluo_GFP_P1.mat';
    file_tif = '190119T132947_Blue_Koehler_Fluo_GFP.tif';
    path_mat = '\\FANTASIA-DS3617\Test_Imaging\2018.11.T1 (Marmoset 132D, Xintrinsic)\M132D-190119-10\';
    load(fullfile(path_mat,file_mat));
    [~,nametemp,~,] = fileparts(file_mat);
    load(fullfile(path_mat,[nametemp(1:end-3),'.mat']));
    I = imread(fullfile(path_mat,file_tif));
    
    case 4 % 96B
    file_mat = '190806T151151_Green_Koehler_Pola_PBS_P1.mat';
    file_tif = '190806T133610_Green_Koehler_Pola_PBS.tif';
    path_mat = '\\FANTASIA-DS3617\Test_Imaging\2018.12.T1 (Marmoset 96B, Xintrinsic)\M96B-190806-13';
    load(fullfile(path_mat,file_mat));
    [~,nametemp,~,] = fileparts(file_mat);
    load(fullfile(path_mat,[nametemp(1:end-3),'.mat']));
    I = imread(fullfile(path_mat,file_tif));
    
    otherwise % browse to open a file
        [file_mat,path_mat] = uigetfile('*.mat','Select a mat file to analyze','\\FANTASIA-DS3617\Test_Imaging\2018.11.T1 (Marmoset 132D, Xintrinsic)');
        load(fullfile(path_mat,file_mat));clc
        [~,nametemp,~,] = fileparts(file_mat);
        load(fullfile(path_mat,[nametemp(1:end-3),'.mat']));
        [file_tif,path_tif] = uigetfile('*.tif','Select a surface image','\\FANTASIA-DS3617\Test_Imaging\2018.11.T1 (Marmoset 132D, Xintrinsic)');
        I = imread(fullfile(path_tif,file_tif));
end


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
if contains(lower(S.SesSoundFile), 'down')
    para.direction = 'down';
else
    para.direction = 'up';
end
clear S P
%% for 80Z, trial-based tonotopy analysis
iStim = [1,4,7];
data_temp       = squeeze(DataMat);
%% View data
opt.ampLimit    = [-0.1, 0.1];
opt.saveON      = 0; 
opt.soundON     = 0;
ViewData_raw(DataMat, para, opt);
%% cycle-based tonotopy analysis

% ===== select frequency component to analysis =====
% period          = S.tm_period; 
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
data_fftamp_mat             = data_fftamp_mat./repmat(data_fftamp_mat(1,:),[L/2+1,1]); % normali zed to mean amplitude
data_fftamp_mean            = mean(data_fftamp_mat,2);
f                           = para.fr*(0:(L/2))/L;

% ========== plot average spectrum ==========
    subplot(2,4,5)
    semilogx(f(2:end),data_fftamp_mean(2:end)), xlim([0 1])
    % label fm and tm periods
%     ind_fm = floor(interp1(f,1:length(f),1/S.fm_period));
    ind = floor(interp1(f,1:length(f),1/period));
    hold on, scatter(f(ind), data_fftamp_mean(ind))
    title('Averaged spectrum')
    % figure, plot(f(2:end),data_fftamp_mat(2:end, sub2ind([para.height, para.width],52, 50))
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
hue     = mod(data_fftagl_one - (delay/20)*2*pi, 2*pi)./(2*pi);

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
saturation_lim  = 0.01; 
saturation      = min(saturation./saturation_lim,1);
map_rgb         = hsv2rgb([hue; saturation; saturation]');
map_rgb         = reshape(map_rgb, [para.height, para.width,3]);
% ========== plot phase map ================
    h2 = subplot(2,4,8);
    imagesc(map_rgb); 
    axis image, colormap(h2, hsv), colorbar
% ===========================================
colorbar('Ticks',0.8.*(para.preStim + para.durStim.*[0:6]./6)./(para.preStim + para.durStim + para.postStim),...
         'TickLabels',{'A4','A5','A6','A7','A8','A9','A10'})

%% plot contours for tonotopy

map_sat         = reshape(saturation, [para.height, para.width]);
mask            = map_sat>0.5;
map_hue = reshape(hue,[para.height,para.width]).*mask;
% figure,imagesc(map_hue), axis image, colormap(hsv)

zindex = 0.8.*(para.durStim.*[0:.5:6]./6)./(para.preStim + para.durStim + para.postStim);
hold all,[M,c] = contour(map_hue, zindex); axis image; c.LineWidth = 2;