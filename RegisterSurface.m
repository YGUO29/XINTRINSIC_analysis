% register surface images across days/sessions
% tonotopy session is the moving image, experiment session if the fixed 
% =========== use: ================== 
% regopt = struct; 
% regopt.mode = 'auto'; % manual or auto
% regopt.manual_method = 'polynomial'; 
% regopt.auto_modality = 'multimodal';
% manually select images: [img_tone_reg, img_exp, tform] = RegisterSurface(para, regopt);
% pre-defined images: [img_tone_reg, img_exp, tform] = RegisterSurface(filename_tone, filename_exp, para, regopt);

function [img_tone_reg, img_exp, tform] = RegisterSurface(varargin)
if length(varargin) == 2 % select images manually
    % open two surface images
    [file,path]     = uigetfile('X:\*.tif','Select a moving image (tonotopy session)');
    filename_tone   = fullfile(path, file);
    [file,path]     = uigetfile('X:\*.tif','Select a fixed image (experiment session)');
    filename_exp    = fullfile(path, file);
    para            = varargin{1};
    opt             = varargin{2};
elseif length(varargin) == 4 % pre-defined input images
    filename_tone   = varargin{1};
    filename_exp    = varargin{2};
    para            = varargin{3};
    opt             = varargin{4};
else
end

% filename_exp    = 'X:\2018.03 T2 (Marmoset 80Z, Xintrinsic, Green & Fluo)\M80Z-180724-13\180724T162507_Blue_Koehler_Fluo_GFP.tif';
% filename_tone   = 'X:\2018.03 T2 (Marmoset 80Z, Xintrinsic, Green & Fluo)\M80Z-180726-14\180726T143800_Blue_Koehler_Fluo_GFP.tif';
% filename_tone   = 'X:\2018.03 T2 (Marmoset 80Z, Xintrinsic, Green & Fluo)\M80Z-180718-14\180718T174201_Blue_Koehler_Fluo_GFP.tif';

% the moving image
img_tone        = imread(filename_tone);
img_tone        = double(imresize(img_tone,para.height./size(img_tone,1)));
img_tone_norm   = (img_tone - min(img_tone(:)))./range(img_tone(:));
% the fixed image
img_exp         = imread(filename_exp);
img_exp         = double(imresize(img_exp,para.height./size(img_exp,1)));
img_exp_norm   = (img_exp - min(img_exp(:)))./range(img_exp(:));

if strcmp(opt.mode, 'manual') % used control points to manully align
    [movingpts, fixedpts] = cpselect(img_tone_norm, img_exp_norm, 'wait', true);
    tform = fitgeotrans(movingpts,fixedpts, opt.manual_method);
else % audomatic alignment
    % registration
    [optimizer,metric]  = imregconfig(opt.auto_modality);
    tform               = imregtform(img_tone, img_exp, 'rigid', optimizer, metric);
end

img_tone_reg        = imwarp(img_tone, tform,'OutputView',imref2d(size(img_exp)));
img_tone_regnorm    = (img_tone_reg - min(img_tone_reg(:)))./range(img_tone_reg(:));

figurex;
subplot(1,2,1), imshowpair(img_exp_norm, img_tone_norm,'Scaling','joint'), title('before registration')
subplot(1,2,2), imshowpair(img_exp_norm, img_tone_regnorm,'Scaling','joint'), title('after registration')

end