% clear all

% load recording session
[file_mat,path_mat] = uigetfile('*.mat','Select a mat file to analyze','X:\');
load(fullfile(path_mat,file_mat));clc
% [~,nametemp,~,] = fileparts(file_mat);
% Sys = load(fullfile(path_mat,[nametemp(1:end-3),'.mat']));
file_parts = strsplit(file_mat,{'_','.'});
ind = find(strcmp(file_parts,'P1'));
nametemp = strjoin(file_parts(1:ind-2),'_');
load(fullfile(path_mat,[nametemp,'.mat']));

% I = double(imresize(I,1/16));
para.nRep =     size(P.ProcDataMat,1);
para.nStim =    size(P.ProcDataMat,2);
para.height =   size(P.ProcDataMat,3);


para.width =    size(P.ProcDataMat,4);
para.nFrame =   size(P.ProcDataMat,5);


para.fr =       P.ProcFrameRate; % frame rate
% para.preStim =  Sys.S.TrlDurPreStim;
% para.durStim =  Sys.S.TrlDurStim;
% para.postStim = Sys.S.TrlDurPostStim;
para.preStim =  S.TrlDurPreStim;
para.durStim =  S.TrlDurStim;
para.postStim = S.TrlDurPostStim;
para.filename = file_mat;
para.pathname = path_mat;
DataMat =       P.ProcDataMat;

% load stimulus information
% if contains(Sys.S.SesSoundFile, 'DMR')
if contains(S.SesSoundFile, 'DMR')
    [file_mat,path_mat] = uigetfile('*.mat','Select a mat file to analyze','D:\=sounds=\DMR');
    load(fullfile(path_mat,file_mat));
    para.tonotopy = 0; % not tonotopy mode
else
    para.tonotopy = 1; % yes, tonotopy mode
end

% if contains(lower(Sys.S.SesSoundFile), 'down')
if contains(lower(S.SesSoundFile), 'down')
    para.direction = 'down';
else
    para.direction = 'up';
end

para.durStim = para.preStim + para.postStim + para.durStim;
para.preStim = 0;
para.postStim = 0;

clear Sys P
size(DataMat)

%% View data
opt.ampLimit    = 0.02.*[-1, 1];
opt.saveON      = 0; 
opt.soundON     = 0;
ViewData_raw(-DataMat, para, opt);
        
%% cycle-based tonotopy analysis
plotON = 1;

% set polarity and delay values
if contains(nametemp,'Fluo')
    para.modality = 'fluo';
    delay = 0;
    data_temp       = squeeze(DataMat);
else
    para.modality = 'intrinsic';
    delay = 3.5;
    data_temp       = -squeeze(DataMat); 
end

% select the frequency component to be analyzed 
if para.tonotopy
    period          = para.preStim + para.durStim + para.postStim;
    rep             = para.nRep;
else
    period          = S.tm_period; 
    rep             = S.tm_cycles;
end

% concatenate repetitions (total frames x total #pixs)
if para.nRep == 1
    data_temp   = permute(data_temp,[3 1 2]);
else
    data_temp   = permute(data_temp,[4 1 2 3]);
end
data_mean       = squeeze(mean(data_temp,1));
nPix            = para.height*para.width;
data_temp       = reshape(data_temp,[para.nRep*para.nFrame, para.height*para.width]);

% average for all pixels
data_temp_mean  = mean(data_temp,2);
% Max             = max(data_temp_mean) + 3*std(data_temp_mean);
% Min             = 2*mean(data_temp_mean) - Max; 
% data_temp_mean  = data_temp(: , sub2ind([para.height, para.width],40,61));

% index for the beginning of each period
ind             = para.fr.*[period: period: floor(period.*rep)];
ind             = [1,ind];
% ===========================================

% =================== FFT =================== 
L               = floor(para.fr * para.nRep * (para.preStim + para.durStim + para.postStim)); %length of signal
T               = 1/para.fr; % sampling period
t               = (0:L-1)*T;
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

freq_comp = floor(interp1(f,1:length(f),1/period));
data_fftagl_mat             = mod( angle(data_fft')', 2*pi );
% data_fftagl_mat             = angle(data_fft')';
% data_fftagl_mat             = angle(data_fft(1:L/2+1,:,:));
% data_fftagl_mat             = data_fftagl_mat + pi.*ones(size(data_fftagl_mat));

% take the corresponding frequency component
data_fftamp_one = squeeze(data_fftamp_mat(freq_comp,:));
data_fftagl_one = squeeze(data_fftagl_mat(freq_comp,:));
% data_fftamp_one = reshape(data_fftamp_one, nPix, 1);
% data_fftagl_one = reshape(data_fftagl_one, nPix, 1);

% Compensate the pseudo-delay
hue     = mod(data_fftagl_one - (delay/period)*2*pi, 2*pi)./(2*pi);
% Reverse the hue for DOWN cycle
if strcmp(para.direction,'down')
    hue = 1 - hue;
end

% cut off a little bit of colorbar, otherwise the two ends are too similar
% and not differentiable
hue_lim = 0.8;
hue = min(hue_lim, hue);
hue_map = reshape(hue,[para.height,para.width]);

saturation      = data_fftamp_one;
saturation_lim  = prctile(saturation,90); 
saturation      = min(saturation./saturation_lim,1);
sat_map         = reshape(saturation, [para.height, para.width]);

map_rgb         = hsv2rgb([hue; saturation; saturation]');
map_rgb         = reshape(map_rgb, [para.height, para.width,3]);

if plotON
    fig = figurex; size_scr = get(0,'ScreenSize'); set(gcf,'position',[1 1 size_scr(3:4)])
    
    subplot(2,4,1:4)
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
    
    % ========== plot average spectrum ==========
    subplot(2,4,5)
    semilogx(f(2:end),data_fftamp_mean(2:end)), xlim([0 1])
    if para.tonotopy
        ind = floor(interp1(f,1:length(f),1/period));
        hold on, scatter(f(ind), data_fftamp_mean(ind))
    else
        % label fm and tm periods
        ind_fm = floor(interp1(f,1:length(f),1/S.fm_period));
        ind_tm = floor(interp1(f,1:length(f),1/S.tm_period));
        hold on, scatter(f([ind_fm, ind_tm]), data_fftamp_mean([ind_fm, ind_tm]))
    end
    title('Averaged spectrum')
    
    % ========== plot phase ====================
    subplot(2,4,6) 
    hist(data_fftagl_one./pi,100);
    title('Phase distribution (uncompensated)')
    
    % ========== plot amplitude map ================
    h1 = subplot(2,4,7);
%     map_amp = reshape(saturation, para.height, para.width); 
    imagesc(sat_map);
    axis image, colormap(h1, gray), colorbar
    
    % ========== plot phase map ================
    h2 = subplot(2,4,8);
    imagesc(map_rgb); 
    axis image, colormap(h2, hsv), colorbar
    
    %     hold on
    %     map_sat         = reshape(saturation, [para.height, para.width]);
    %     mask            = map_sat>0.4;
    %     map_hue         = reshape(hue,[para.height,para.width]).*mask;
    %     zindex = (2.7+12*0.2: 12*0.2: 20-2.7)./(para.preStim + para.durStim + para.postStim);
    %     hold all,[M,c] = contour(map_hue, zindex); axis image; c.LineWidth = 3;     
    %     colormap(h2, hsv), caxis([2.7 20-2.7]./20)
    % ===========================================
end

% % ============== Sine fm =============
%      colorbar('Ticks',0.8.*(para.preStim + para.durStim.*[0:8]./8)./(para.preStim + para.durStim + para.postStim),...
%          'TickLabels',{'4','6.8','8','6.8','4','1.28','0.125','1.28','4'})
% ============== SineLog fm =============
% colorbar('Ticks',0.8.*(para.preStim + para.durStim.*[0:8]./8)./(para.preStim + para.durStim + para.postStim),...
%          'TickLabels',{'1','4.35','8','4.35','1','0.23','0.125','0.23','1'})
 % ============== Sine tm =============
% colorbar('Ticks',0.8.*(para.preStim + para.durStim.*[0:8]./8)./(para.preStim + para.durStim + para.postStim),...
%          'TickLabels',{'0','22.6','32','22.6','0','-22.6','-32','-22.6','0'})

%% another way for contour (non-overlap)

map_sat         = reshape(saturation, [para.height, para.width]);
mask            = map_sat>0.4;
map_hue         = reshape(hue,[para.height,para.width]).*mask;
zindex = (2.7: 12*0.2: 20-2.7)./(para.preStim + para.durStim + para.postStim);

quant_hue       = imquantize(map_hue,zindex);
figurex;
bd = cell(1);
for i = 2:max(quant_hue(:))
    temp = quant_hue==i;
%     bd{i-1} = boundarymask(temp);
    contour(temp)
end
%% plot two amplitude maps together
map_temp = zeros(para.height, para.width, 3);
saturation_lim = 0.03;
period          = S.fm_period; % select frequency component to analysis
freq_comp = floor(interp1(f,1:length(f),1/period));
saturation = squeeze(data_fftamp_mat(freq_comp,:));
saturation = min(saturation./saturation_lim,1);
saturation = reshape(saturation, para.height, para.width); 
map_temp(:,:,1) = saturation;

saturation_lim = 0.02;
period          = S.tm_period; % select frequency component to analysis
freq_comp = floor(interp1(f,1:length(f),1/period));
saturation = squeeze(data_fftamp_mat(freq_comp,:));
saturation = min(saturation./saturation_lim,1);
saturation = reshape(saturation, para.height, para.width); 

map_temp(:,:,2) = saturation;

figure, imagesc(map_temp), 
axis image
     

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