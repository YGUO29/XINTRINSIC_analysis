function [newT] = getRegisteredTonotopy(T, para, animal, mode)

switch animal
    case '102D'
        filename_exp = 'U:\2019.09.T1 (Marmoset 102D, Xintrinsic)\M102D-201226-11\201226T114949_Blue_Koehler_Fluo_GFP.tif';
    %     filename_exp = 'U:\2019.09.T1 (Marmoset 102D, Xintrinsic)\M102D-210205-09\210205T095705_Blue_Koehler_Fluo_GFP.tif'; % nat 168
        filename_exp = 'U:\2019.09.T1 (Marmoset 102D, Xintrinsic)\M102D-210121-09\210121T122526_Blue_Koehler_Fluo_GFP.tif'; % voc major
    %     filename_exp = 'U:\2019.09.T1 (Marmoset 102D, Xintrinsic)\M102D-210210-09\210210T100434_Blue_Koehler_Fluo_GFP.tif';
    case '80Z'
        filename_exp    = 'Y:\2018.03_M80Z\M80Z-180724-13\180724T162507_Blue_Koehler_Fluo_GFP.tif';
        % filename_tone   = 'X:\2018.03 T2 (Marmoset 80Z, Xintrinsic, Green & Fluo)\M80Z-180726-14\180726T143800_Blue_Koehler_Fluo_GFP.tif';
        filename_tone   = 'Y:\2018.03_M80Z\M80Z-180726-14\180726T173601_Blue_Koehler_Fluo_GFP.tif';
        % intrinsic tonotopy session: 180621T142313_Blue_Koehler_Fluo_GFP.rec
        % calcium tonotopy session: 180726T144709_Blue_Koehler_Fluo_GFP.rec
    case '132D'
    %     filename_exp    = 'X:\2018.11.T1 (Marmoset 132D, Xintrinsic)\M132D-190224-09\190224T095029_Blue_Koehler_Fluo_GFP.tif';
        filename_exp    = 'Y:\2018.11_M132D\M132D-190119-10\190119T132947_Blue_Koehler_Fluo_GFP.tif';
        filename_tone   = 'Y:\2018.11_M132D\M132D-190409-10\190409T133232_Blue_Koehler_Fluo_GFP.tif';
    %     filename_tone   = 'X:\2018.11.T1 (Marmoset 132D, Xintrinsic)\M132D-190326-12\190326T142236_Blue_Koehler_Fluo_GFP.tif';
    %     % tonotopy session: 190326T134033_Blue_Koehler_Fluo_GFP.rec
    case '126D'
    %% 126D (same session)
        filename_exp    = 'X:\2019.05.T1 (Marmoset 126D, Xintrinsic)\M126D-191023-08\191023T133842_Green_Koehler_Pola_PBS.tif';
        filename_tone   = 'X:\2019.05.T1 (Marmoset 126D, Xintrinsic)\M126D-191023-08\191023T133842_Green_Koehler_Pola_PBS.tif';
    otherwise
end

%% register two surface image
regopt = struct; 
regopt.mode = mode; % manual or auto
regopt.manual_method = 'pwl'; 
regopt.auto_method = 'similarity'; % translation / rigid / similarity / affine
regopt.auto_modality = 'multimodal';
[img_tone_reg, img_exp, tform] = RegisterSurface(filename_tone, filename_exp, para, regopt);
% [img_tone_reg, img_exp, tform] = RegisterSurface(para, regopt);

%% step3 (for cycle-based tonotopy): move tonotopy to match experiment session (with tform)
% saturation_reg = imwarp(sat_map, tform,'OutputView',imref2d(size(img_exp)));
% hue_reg = imwarp(hue_map, tform,'OutputView',imref2d(size(img_exp)));
% figurex;
% subplot(1,2,1), imshowpair(sat_map, saturation_reg,'Scaling','joint')
% subplot(1,2,2), imshowpair(hue_map, hue_reg,'Scaling','joint')
% map_rgb_reg = imwarp(map_rgb, tform,'OutputView',imref2d(size(img_exp)));

%% step3 (for trial-based tonotopy): move tonotopy to match experiment session (with tform) 
newT = T;
for i = 1:3
newT.tonotopy{i} = imwarp(T.tonotopy{i}, tform, 'OutputView', imref2d(size(img_exp)));
end
newT.sat_map = imwarp(T.sat_map, tform, 'OutputView', imref2d(size(img_exp)));
T.hue_map(isnan(T.hue_map)) = 0; % a few values might be NAN 
newT.hue_map = imwarp(T.hue_map, tform, 'OutputView', imref2d(size(img_exp)));

figurex;
subplot(1,2,1), imagesc(T.tonotopy{3}), axis image, axis off
title('original')
subplot(1,2,2), imagesc(newT.tonotopy{3}), axis image, axis off
title('transformed')

