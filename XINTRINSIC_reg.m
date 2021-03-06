% registration with matlab function
function [DataMat_reg, Yreg] = XINTRINSIC_reg(Y, para, opt)
addpath(genpath('D:\SynologyDrive\=code=\FANTASIA-NoRMCorre'));

%% define template
% bnd:          number of pixels to be exluded to avoid NaN effects 
%                   [x_beg,x_end,y_be,y_end,z_beg,z_end]
bnd = 5;
[cY,mY,~] = motion_metrics(Y, bnd);
% cY:           correlation coefficient of each frame with the mean
% mY:           mean image
% ng:           norm of gradient of mean image

% fixed = Y(:,:,1);
fixed = mY;
nRegFrames = length(opt.ind_reg);
% movingRegistered = moving;
[optimizer,metric] = imregconfig('multimodal'); % monomodal or multimodal

Yreg = Y;
f = waitbar(0, 'registering...');
for i = 1:nRegFrames
% for i = 280:340

tstart = tic;
    waitbar(i/nRegFrames, f, ['registering ', num2str(i), '/', num2str(nRegFrames)])
    tform = imregtform(Y(:,:,opt.ind_reg(i)),fixed,'rigid', optimizer, metric);
    Yreg(:,:,opt.ind_reg(i)) = imwarp(Y(:,:,opt.ind_reg(i)),tform,'OutputView',imref2d(size(fixed)));
end
time_perframe = toc(tstart)/300
close(f)

[cY1,mY1,~] = motion_metrics(Yreg,bnd);
if opt.plotON
    % plot correlation coefficients
    figurex;
    subplot(2,2,[1,2])
    t = 1:length(cY1);
    plot(t, cY), xlim([0, max(t)]), ylim([min(cY), 1]), hold on
    plot(t, cY1), xlim([0, max(t)]), ylim([min(cY1), 1])
    
    subplot(2,2,3)
    imshow(mY,[]), title('before')
    
    subplot(2,2,4)
    imshow(mY1,[]), title('after')
end

DataMat_reg = reshape(Yreg, para.height, para.width, para.nFrame, para.nStim, para.nRep);
[~,ind] = sort(para.order);
tic, DataMat_reg = DataMat_reg(:,:,:,ind,:); tReorder = toc % re-arrange according to the experiment order
tic, DataMat_reg = permute(DataMat_reg,[5, 4, 1, 2, 3]); tPermute = toc % DataMat_reg = [height, width, frames, trial, rep]

end