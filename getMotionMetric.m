function [output, opt] = getMotionMetric(DataMat, para, opt)
addpath(genpath('D:\SynologyDrive\=code=\FANTASIA-NoRMCorre'));

if length(size(DataMat)) == 5
    % put all frames in a 3d matrix: Y
    frames = permute(DataMat,[3,4,5,2,1]);% DataMat = [rep, trial, height, width, frams]
    % re-arrange according to the experiment order, so that the time points are continuous
%     frames = frames(:,:,:, para.order, opt.reps); % frames = [height,
%     width, frames, trial, rep] %(ERROR: para.order does not change for each repetition, not necessary)
    frames = frames(:,:,:, :, opt.reps); % frames = [height, width, frames, trial, rep] 
    Y = reshape(frames, para.height, para.width, para.nFrame*para.nStim*length(opt.reps)); 
else % for registering across sessions
    Y = DataMat;
end

% test a few frames
% frames = DataMat(:, [1 80 180], :,:,:);
% frames = permute(frames,[3,4,1,2,5]);
% Y = reshape(frames, [para.height, para.width, size(frames,3)*size(frames,4)*size(frames,5)]); 

% visualize frames
% figure,
% for i = 1:2000
%     imshow(Y(:,:,i),[]),title(num2str(i))
%     pause(0.1)
% end

% calculate metrics 
Y = single(Y);
Y = Y - min(Y(:));

% bnd:          number of pixels to be exluded to avoid NaN effects 
%                   [x_beg,x_end,y_be,y_end,z_beg,z_end]

% cY:           correlation coefficient of each frame with the mean
% mY:           mean image
% ng:           norm of gradient of mean image
[cY,mY,~] = motion_metrics(Y, opt.bnd); % cY: correlation coefficients between frames

%% define a threshold for "moving frames"
opt.cY_threshold = prctile(cY, opt.th_prctile);
if opt.plotON
    % plot correlation coefficients
    figurex;
%     t = 1:size(frames,3)*size(frames,4)*size(frames,5);
    t = 1:size(Y, 3);
    plot(t, cY), xlim([0, max(t)]), ylim([min(cY), 1]), title('Original')
    % plot histogram and motion correction threshold of correlation coefficients
    figurex; 
    hist(cY,100);
    [counts, ~] = hist(cY,100);
    hold on
    line([opt.cY_threshold, opt.cY_threshold], [0 max(counts)])
    xlabel('cY'), ylabel('counts')
end

opt.ind_reg = find(cY < opt.cY_threshold);

%%
output.Y = Y;
output.cY = cY;
output.mY = mY;
end