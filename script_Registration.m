% registration across sessions / correct motion artifacts
%% check alignment for all repetitions
opt_reg.reps = 4:para.nRep;
opt_reg.plotON = 1;
opt_reg.th_prctile = 30; % define x percentile of total frames will not be registered (x% will be registered)
opt_reg.bnd = 5; % boundary pixels excluded from motion metrics measurement

% check DataMat (5d matrix)
[output, opt_reg] = getMotionMetric(DataMat, para, opt_reg); 

% check Y (3d matrix)
% [output, opt_reg] = getMotionMetric(reshape(Yreg_new, size(Yreg_new,1), size(Yreg_new,2), size(Yreg_new,3)*size(Yreg_new,4))...
%     , para, opt_reg); 
%% method 1 - register the entire session all together 
%% perform registration with NoRMCorre
opt_reg.plotON = 1;
opt_reg.mode = 'rigid';
opt_reg.mask = 0;
opt_reg.reps = 1:para.nRep;
opt_reg.regreps = 0; %registration across repetitions (register mean images obtained from repetitions)
% opt_reg.template = mY; 
tstart = tic;
[output_reg, opt_reg] = XINTRINSIC_reg_NoRMCorre(output.Y, para, opt_reg);
time.registration = toc(tstart);
DataMat_reg = output_reg.DataMat_temp;
[output, opt_reg] = getMotionMetric(DataMat_reg, para, opt_reg);


%% method 2 - register session by session (usually recorded on different days)
%% get motion metrics for each session
% initialize parameters
opt_reg = struct;
DataMat_reg = DataMat;
Y = zeros(para.height, para.width, para.nStim * para.nFrame, para.nRep);% frames (H, W, trials x frames, reps)
Yreg = Y;
cY = zeros(para.nStim * para.nFrame, para.nRep);% corrcoef (trials x frames, reps)
cYreg = cY;
mY = zeros(para.height, para.width, para.nRep);% mean image (H, W, reps)

% register within each repetition
for i = 1:para.nRep
    % get motion metrics, determine which frames to register (index stored
    % in opt_reg.ind_reg)
    opt_reg.reps = i;
    opt_reg.plotON = 0;
    opt_reg.th_prctile = 30; % define x percentile of total frames will not be registered (x% will be registered)
    opt_reg.bnd = 5; % boundary pixels excluded from motion metrics measurement
    opt_reg.regreps = 0;
    tstart = tic;
    [output, opt_reg] = getMotionMetric(DataMat, para, opt_reg); 
    % (updated opt_reg.cY_threshold and opt_reg.ind_reg)
    Y(:,:,:,i) = output.Y;
    cY(:,i) = output.cY;   
    time.motionmetric = toc(tstart);
    
    % registration (selected frames)
    opt_reg.mode = 'rigid';
    opt_reg.mask = 0;
    % opt_reg.template = mY;
    tstart = tic;
    [output, opt_reg] = XINTRINSIC_reg_NoRMCorre(output.Y, para, opt_reg);
    DataMat_reg(i, :, :, :, :) = output.DataMat_temp; % (H, W, frame, stim, 1)
    Yreg(:,:,:,i) = output.Yreg;
    cYreg(:,i) = output.cYreg;
    mY(:,:,i) = output.mYreg;
    time.registration = toc(tstart);
end

% temp
% transform DataMat to Yf
frames = permute(DataMat_reg, [3,4,5,2,1]);% DataMat = [rep, trial, height, width, frams]
% re-arrange according to the experiment order, so that the time points are continuous
Yreg = reshape(frames, size(frames, 1), size(frames, 2), size(frames, 3)*size(frames, 4), size(frames, 5)); 

for i = 1:size(frames, 4)
    mY(:,:,i) = squeeze(mean(frames(:,:,:, i), 3));
end
%% register mY across reps, get shifts
opt_reg.ind_reg = 1:size(mY, 3);
opt_reg.plotON = 1;
opt_reg.th_prctile = 100;
opt_reg.regreps = 1; %registration across repetitions (register mean images obtained from repetitions)
opt_reg.template = mY(:,:,1);
opt_reg.mask = 0;
opt_reg.mode = 'rigid';
opt_reg.bnd = 5; % boundary pixels excluded from motion metrics measurement

[output, opt_reg] = XINTRINSIC_reg_NoRMCorre(mY, para, opt_reg);
shifts = output.shifts;
options = output.options;
mYreg = output.Yreg;


% apply shifts to all frames for each rep
Yreg_new = Yreg;
for i = 2:para.nRep
    shifts_temp = shifts(i);
    shifts_group(1:size(Yreg, 3)) = shifts_temp;
    tstart = tic;
    Yreg_new(:,:,:,i) = apply_shifts(squeeze(Yreg(:,:,:,i)), shifts_group, options);
%     Yreg_new = apply_shifts('Yreg_test.mat', shifts(2), options);
    time.applyshift = toc(tstart);
end
%% visualize 2 sets of frames
% Y_ori = squeeze(Y(:,:,:,12)); 
% Y_reg = squeeze(Yreg(:,:,:,12)); 
% Y_ori = mY;
% Y_reg = mYreg;

nnY = quantile(Y_ori(:),0.005);
mmY = quantile(Y_ori(:),0.995);
T = size(Y_ori,3);
figurex([1250         351        1848         545]);
for t = 101000:1020000
    subplot(121);imagesc(Y_ori(:,:,t),[nnY,mmY]); xlabel('raw data','fontsize',14,'fontweight','bold'); axis equal; axis tight;
    title(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14); colormap('bone')
    subplot(122);imagesc(Y_reg(:,:,t),[nnY,mmY]); xlabel('non-rigid corrected','fontsize',14,'fontweight','bold'); axis equal; axis tight;
    title(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14); colormap('bone')
    set(gca,'XTick',[],'YTick',[]);
    drawnow;
    pause(0.001);
end
%% perform registration with MATLAB functions
% [DataMat_reg, Yreg] = XINTRINSIC_reg(Y, para, opt_reg);
% [~, cYreg, ~] = getMotionMetric(DataMat_reg, para, opt_reg);
