clear opt_read
opt_read.animal = '132D'; 
opt_read.session = 'Tone';
opt_read.modal = 'Calcium';
opt_read.date = '190409';
% opt_read.subset = '50dB'; % only read subset of the data (comment this line if reading all data)
[DataMat, para] = readSession(opt_read);

DataMat = double(DataMat);
opt = struct;
[X, DataMat_norm, mov_rel_sep] = getX(DataMat, para, opt);

% ============== get tonotopy ==================
para = getTonotopyParameters(para, '80Z');
T = getTonotopy(X, para);

% register if necessary
T = getRegisteredTonotopy(T, para, '132D', 'manual'); 
% for 132D: use manual alignment

% get tonotopy contour
M.saturation = T.sat_map; % used for masking
M.hue = T.hue_map; % used for plotting contour
% chose a map to overlay the contour on
M.map_rgb = T.tonotopy{3};
para.sat_lim = 0.02;
para.zindex = T.zindex;
ct = getContour(M, para);

figurex; 
imagesc(M.map_rgb); 
% imagesc(map_hue); 
axis image, axis off
plotContour(ct)

clear opt_read
opt_read.animal = '132D'; 
opt_read.session = 'Nat';
opt_read.modal = 'Calcium';
opt_read.date = '190119';
% opt_read.subset = 'Nat'; % only read subset of the data (comment this line if reading all data)
[DataMat, para] = readSession(opt_read);

DataMat = double(DataMat);
opt = struct;
[X, DataMat_norm, mov_rel_sep] = getX(DataMat, para, opt);

% ============= pixel-wise analysis =====================
load('D:\SynologyDrive\=data=\F_halfcosine_marm_NatJM.mat') 
clear output
opt_corr.CV = 0;
opt_corr.comp = 0;
opt_corr.regress_type = 'regular';
output = getCorrcoef(X, F, opt_corr);

% plot coefficients
para.mirror = 0;
para.ct = ct;
plotRegressCoef(output, F, para)
plotRegressCoef2(output, F, para)
plotRegressRsquared(output, F, para)
plotRegressOverlay(output, para)

plotRegressCorr(output, F, T, para, '80Z')