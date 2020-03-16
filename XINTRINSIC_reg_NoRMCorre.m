% modified by Yueqi 2018/8
function XINTRINSIC_reg_NoRMCorre(DataMat, I, para, opt)
addpath(genpath('D:\=code=\FANTASIA-NoRMCorre'));

tic, DataMat = permute(DataMat,[3,4,5,2,1]); time.permute = toc % DataMat = [rep, trial, height, width, frams]
tic, DataMat = DataMat(:,:,:,para.order,:); time.reorder = toc % re-arrange according to the experiment order
Y = reshape(DataMat, para.height, para.width, para.nFrame*para.nStim*para.nRep); 

if opt.mask
    figure, imshow(I,[])
    % h = images.roi.Circle(gca,'Center',[floor(para.width/2) floor(para.height/2)],'Radius',floor(para.height/2)); 
    h = images.roi.Rectangle(gca,'Position',[floor(para.width/2)-floor(para.height/2), 1,...
                                             floor(para.height), floor(para.height)]); % xmin, ymin, width, height
    % h = images.roi.Ellipse(gca,'Center',[floor(para.width/2) floor(para.height/2)],'Semiaxes',[40 20]); 
    % h = images.roi.Polygon(gca,'Position',[1 1; 1 para.height - 30; 30 para.height; para.width para.height; para.width 1]); 
    title('Press Enter after the position is adjusted')
    pause
    Xmin = round(h.Position(1)); 
    Ymin = round(h.Position(2));
    Xmax = Xmin + round(h.Position(3));
    Ymax = Ymin + round(h.Position(4));
    Zmin = 1;
    Zmax = size(Y,3);
else
    Xmin = 1; Ymin = 1;
    Xmax = para.width;
    Ymax = para.height;
    Zmin = 1;
    Zmax = size(Y,3);
end

Y = single(Y);                 % convert to single precision 
Y = Y - min(Y(:));

%% set parameters (first try out rigid motion correction)
template = Y(:,:,1); % use the first frme as template
options_rigid = NoRMCorreSetParms(...
    'd1',           size(Y,1),...
    'd2',           size(Y,2),...
    'bin_width',    10,...
    'max_shift',    10,...
    'us_fac',       20,...
    'init_batch',   200,...
    'correct_bidir',0); % do not perform bi-directional scanning correction (avoid zigzag artifact)

% perform motion correction
tic; [Yreg_part,shifts1,template1,options_rigid] = normcorre(Y(Ymin:Ymax, Xmin:Xmax, Zmin:Zmax),options_rigid, template); tRegistration = toc
if opt.mask
    Yreg = apply_shifts(Y,shifts1,options_rigid);
else
    Yreg = Yreg_part;
end

DataMat_reg = reshape(Yreg, para.height, para.width, para.nFrame, para.nStim, para.nRep);
[~,ind] = sort(para.order);
tic, DataMat_reg = DataMat_reg(:,:,:,ind,:); tReorder = toc % re-arrange according to the experiment order
tic, DataMat_reg = permute(DataMat_reg,[5, 4, 1, 2, 3]); tPermute = toc % DataMat_reg = [height, width, frames, trial, rep]

    
% save file as tiff
% M1 = uint16(M1);
% saveastiff(M1,'D:\=data=\80Z_imaging\img_2p\CaImAn_Results\20180503T121108_nobidir.tif');
% saveastiff(template1,'D:\=data=\80Z_imaging\img_2p\CaImAn_Results\20180503T121108_template.tif');
%% now try non-rigid motion correction (also in parallel)
% clear large data in memory, to avoid out of memory problems
options_nonrigid = NoRMCorreSetParms(...
    'd1',size(Y,1),...
    'd2',size(Y,2),...
    'grid_size',[16, 16],...
    'mot_uf',4,...
    'bin_width',200,...
    'max_shift',10,...
    'max_dev',3,...
    'us_fac',50,...
    'init_batch',200,...
    'correct_bidir',0); % do not perform bi-directional scanning correction (avoid zigzag artifact)

% perform motion correction
tic; [M2,shifts2,template2,options_nonrigid] = normcorre_batch(Y,options_nonrigid); toc
% M2 = uint16(M2);
% saveastiff(M2,'D:\=data=\80Z_imaging\img_2p\CaImAn_Results\20180503T121108_nonrigid_nonbidir.tif');
% saveastiff(template2,'D:\=data=\80Z_imaging\img_2p\CaImAn_Results\20180503T121108_template2.tif');
%% compute metrics

nnY = quantile(Y(:),0.005);
mmY = quantile(Y(:),0.995);

[cY,mY,vY] = motion_metrics(Y,10);
[cM1,mM1,vM1] = motion_metrics(Yreg,10);
% [cM2,mM2,vM2] = motion_metrics(M2,10);
T = length(cY);
%% plot metrics
figure;

    ax1 = subplot(2,2,1); imagesc(mY,[nnY,mmY]);  axis equal; axis tight; axis off; title('mean raw data','fontsize',14,'fontweight','bold')
    ax2 = subplot(2,2,2); imagesc(mM1,[nnY,mmY]);  axis equal; axis tight; axis off; title('mean rigid corrected','fontsize',14,'fontweight','bold')
%     ax3 = subplot(2,3,3); imagesc(mM2,[nnY,mmY]); axis equal; axis tight; axis off; title('mean non-rigid corrected','fontsize',14,'fontweight','bold')
    subplot(2,2,3); plot(1:T,cY,1:T,cM1); legend('raw data','rigid'); title('correlation coefficients','fontsize',14,'fontweight','bold')
    subplot(2,2,4); scatter(cY,cM1); hold on; plot([0.9*min(cY),1.05*max(cM1)],[0.9*min(cY),1.05*max(cM1)],'--r'); axis square;
        xlabel('raw data','fontsize',14,'fontweight','bold'); ylabel('rigid corrected','fontsize',14,'fontweight','bold');
%     subplot(2,3,6); scatter(cM1,cM2); hold on; plot([0.9*min(cY),1.05*max(cM1)],[0.9*min(cY),1.05*max(cM1)],'--r'); axis square;
%         xlabel('rigid corrected','fontsize',14,'fontweight','bold'); ylabel('non-rigid corrected','fontsize',14,'fontweight','bold');
    linkaxes([ax1,ax2],'xy')
%% plot shifts        

shifts_r = squeeze(cat(3,shifts1(:).shifts));
% shifts_nr = cat(ndims(shifts2(1).shifts)+1,shifts2(:).shifts);
% shifts_nr = reshape(shifts_nr,[],ndims(Y)-1,T);
% shifts_x = squeeze(shifts_nr(:,1,:))';
% shifts_y = squeeze(shifts_nr(:,2,:))';

% patch_id = 1:size(shifts_x,2);
% str = strtrim(cellstr(int2str(patch_id.')));
% str = cellfun(@(x) ['patch # ',x],str,'un',0);

figure;
    ax1 = subplot(311); plot(1:T,cY,1:T,cM1); legend('raw data','rigid'); title('correlation coefficients','fontsize',14,'fontweight','bold')
            set(gca,'Xtick',[])
    ax2 = subplot(312); plot(shifts_r(:,1),'k','linewidth',2); title('displacements along x','fontsize',14,'fontweight','bold')
            set(gca,'Xtick',[])
    ax3 = subplot(313); plot(shifts_r(:,2),'k','linewidth',2); title('displacements along y','fontsize',14,'fontweight','bold')
            xlabel('timestep','fontsize',14,'fontweight','bold')
    linkaxes([ax1,ax2,ax3],'x')

%% plot a movie with the results

figure;
for t = 1336:2924
    subplot(121);imagesc(Y(:,:,t),[nnY,mmY]); xlabel('raw data','fontsize',14,'fontweight','bold'); axis equal; axis tight;
    title(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14); colormap('bone')
    subplot(122);imagesc(Yreg(:,:,t),[nnY,mmY]); xlabel('non-rigid corrected','fontsize',14,'fontweight','bold'); axis equal; axis tight;
    title(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14); colormap('bone')
    set(gca,'XTick',[],'YTick',[]);
    drawnow;
    pause(0.01);
end

end