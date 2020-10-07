% registration across sessions / correct motion artifacts

%% check alignment across sessions
opt_reg.plotON = 1;
opt_reg.th_prctile = 30; % define x percentile of total frames will not be registered (x% will be registered)
opt_reg.bnd = 5; % boundary pixels excluded from motion metrics measurement
tstart = tic;
[Y, cY, opt_reg] = getMotionMetric(DataMat, para, opt_reg);
t_registration = toc(tstart);
%% perform registration with MATLAB functions
% [DataMat_reg, Yreg] = XINTRINSIC_reg(Y, para, opt_reg);
% [~, cYreg, ~] = getMotionMetric(DataMat_reg, para, opt_reg);
%% perform registration with NoRMCorre
opt_reg.mode = 'nonrigid';
opt_reg.mask = 0;
[DataMat_reg, Yreg] = XINTRINSIC_reg_NoRMCorre(Y, para, opt_reg);
[~, cYreg, ~] = getMotionMetric(DataMat_reg, para, opt_reg);
