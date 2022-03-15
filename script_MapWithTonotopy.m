% plot a map from an experiment together with tonotopy
%% step1: process a tonotopy session
% go to script_CycleBased to get sat_map 
% or use trial-based tonotopy (script_TonotopyTrialBased.m)

%% step2: load surface images from tonotopy session and experiment session, and register to get "tform"

    %% 102D,
    filename_tone = 'U:\2019.09.T1 (Marmoset 102D, Xintrinsic)\M102D-201226-11\201226T114949_Blue_Koehler_Fluo_GFP.tif';
%     filename_exp = 'U:\2019.09.T1 (Marmoset 102D, Xintrinsic)\M102D-210205-09\210205T095705_Blue_Koehler_Fluo_GFP.tif'; % nat 168
    filename_exp = 'U:\2019.09.T1 (Marmoset 102D, Xintrinsic)\M102D-210121-09\210121T122526_Blue_Koehler_Fluo_GFP.tif'; % voc major
%     filename_exp = 'U:\2019.09.T1 (Marmoset 102D, Xintrinsic)\M102D-210210-09\210210T100434_Blue_Koehler_Fluo_GFP.tif';
    %% 80Z, tonotopy and pure tone (trial based)
    filename_exp    = 'Y:\2018.03_M80Z\M80Z-180724-13\180724T162507_Blue_Koehler_Fluo_GFP.tif';
    % filename_tone   = 'X:\2018.03 T2 (Marmoset 80Z, Xintrinsic, Green & Fluo)\M80Z-180726-14\180726T143800_Blue_Koehler_Fluo_GFP.tif';
    filename_tone   = 'Y:\2018.03_M80Z\M80Z-180726-14\180726T173601_Blue_Koehler_Fluo_GFP.tif';
    % intrinsic tonotopy session: 180621T142313_Blue_Koehler_Fluo_GFP.rec
    % calcium tonotopy session: 180726T144709_Blue_Koehler_Fluo_GFP.rec
    %% 132D
%     filename_exp    = 'X:\2018.11.T1 (Marmoset 132D, Xintrinsic)\M132D-190224-09\190224T095029_Blue_Koehler_Fluo_GFP.tif';
    filename_exp    = 'U:\2018.11.T1 (Marmoset 132D, Xintrinsic)\M132D-190119-10\190119T132947_Blue_Koehler_Fluo_GFP.tif';
    filename_tone   = 'U:\2018.11.T1 (Marmoset 132D, Xintrinsic)\M132D-190409-10\190409T133232_Blue_Koehler_Fluo_GFP.tif';
%     filename_tone   = 'X:\2018.11.T1 (Marmoset 132D, Xintrinsic)\M132D-190326-12\190326T142236_Blue_Koehler_Fluo_GFP.tif';
%     % tonotopy session: 190326T134033_Blue_Koehler_Fluo_GFP.rec

    %% 126D (same session)
    filename_exp    = 'X:\2019.05.T1 (Marmoset 126D, Xintrinsic)\M126D-191023-08\191023T133842_Green_Koehler_Pola_PBS.tif';
    filename_tone   = 'X:\2019.05.T1 (Marmoset 126D, Xintrinsic)\M126D-191023-08\191023T133842_Green_Koehler_Pola_PBS.tif';
    
    %% register two surface image
    regopt = struct; 
    regopt.mode = 'auto'; % manual or auto
    regopt.manual_method = 'pwl'; 
    regopt.auto_method = 'similarity'; % translation / rigid / similarity / affine
    regopt.auto_modality = 'multimodal';
    [img_tone_reg, img_exp, tform] = RegisterSurface(filename_tone, filename_exp, para, regopt);
    % [img_tone_reg, img_exp, tform] = RegisterSurface(para, regopt);

%% step3 (for cycle-based tonotopy): move tonotopy to match experiment session (with tform)
saturation_reg = imwarp(sat_map, tform,'OutputView',imref2d(size(img_exp)));
hue_reg = imwarp(hue_map, tform,'OutputView',imref2d(size(img_exp)));
figurex;
subplot(1,2,1), imshowpair(sat_map, saturation_reg,'Scaling','joint')
subplot(1,2,2), imshowpair(hue_map, hue_reg,'Scaling','joint')
map_rgb_reg = imwarp(map_rgb, tform,'OutputView',imref2d(size(img_exp)));
%% step3 (for trial-based tonotopy): move tonotopy to match experiment session (with tform) 
map_rgb_reg = imwarp(tonotopy{3}, tform, 'OutputView', imref2d(size(img_exp)));
figurex;
subplot(1,2,1), imagesc(tonotopy{3}), axis image, axis off
title('original')
subplot(1,2,2), imagesc(map_rgb_reg), axis image, axis off
title('transformed')
%% step3 (for cycle-based): if no registration needed
% for cycle based tonotopy:
saturation_reg = sat_map;
hue_reg = hue_map;
map_rgb_reg = map_rgb;

%% step3 (for trial-based): for trial-based tonotopy:
saturation_reg = T.sat_map;
hue_reg = T.hue_map;
map_rgb_reg = T.tonotopy{3};

%% step 4.1: get contour for tonotopy (plot contours with tonotopy to double check)
% stimulus parameter for sound "TonePipSeq_A4-A10_73x0.2s(0.2s)@1.0st_(2.7+14.6+2.7)s"
M.saturation = saturation_reg; % used for masking
M.hue = hue_reg; % used for plotting contour
% chose a map to overlay the contour on
M.map_rgb = map_rgb_reg;
para.sat_lim = 0.02;
para.zindex = T.zindex;
ct = getContour(M, para);
%% step 4.2: get contour for components
% see function: plotTonotopyWithComponents
% input: a comp*pixel matrix (weight matrix)
close all
iK = 1;
nComp = size(Ws{iK},1);
ct = getContour_compmponent(Ws{iK}, comp{iK});

% verification: plot components contour with components
figurex;
for iComp = 1:nComp
    ax = subplot(1, nComp, iComp);
    max_amp = max(max( abs(comp{iK}(:,:,iComp)) ));
    imagesc(comp{iK}(:,:,iComp), [-max_amp, max_amp]), axis image, axis off, colormap(jet)
    if para.mirror
        set(gca,'XDir','reverse')
    else 
    end
    plotContour(ct(iComp),'RdBu');
    
end

% align: plot components contour with tonotopy
figurex;
for iComp = 1:nComp
    subplot(1, nComp, iComp)
    imagesc(T.tonotopy{3}); axis image, axis off
    if para.mirror
        set(gca,'XDir','reverse')
    else 
    end
    plotContour(ct(iComp), 'RdBu');
end
%% step5: prosess the experimental session
% go to script_TrialBased or script_Cyclebased
% get comp
%% step6: plot experimental session with tonotopy contours
figurex;
imagesc(map_rgb_reg); axis image
% colormap(jet)
plotContour(ct);

%% another way to plot tonotopy & components together - use components as a transparency mask
comp_mask = comp;
for k = 1:K
%     comp_mask{k} = -comp{k};
    comp_mask{k}(comp_mask{k} < 0) = 0;
    comp_mask{k} = comp_mask{k}./max(comp_mask{k}(:));
    comp_mask{k} = 1 - comp_mask{k}; % convert to transparency (0 = completely transparent)
end

%% overlay component mask on tonotopy
p = [1, K];
figurex;
set(gcf,'color','w','position', [1323, 380, 1365, 192]);
ha = tight_subplot(p(1),p(2),[.01 .01],[.1 .01],[.01 .01]); 
for k = 1:K
    axes(ha(k));
    % imagesc(zeros(size(comp_mask{3})), 'AlphaData', comp_mask{3}), colormap('gray')
    imagesc(map_rgb_reg); 
    hold on
    h = imagesc(zeros(size(comp_mask{k})), [0 1]); colormap(gray)
    set(h, 'AlphaData', comp_mask{k});
    axis image
    axis off
end