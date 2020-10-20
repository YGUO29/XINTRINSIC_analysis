% modified by Yueqi 2018/9
% function XINTRINSIC_reg_NoRMCorre(DataMat, I, para, opt)
% modified by Yueqi 2020/9
function [output, opt] = XINTRINSIC_reg_NoRMCorre(Y, para, opt)
addpath(genpath('D:\SynologyDrive\=code=\FANTASIA-NoRMCorre'));


[cY, mY,~] = motion_metrics(Y, opt.bnd);
if isfield(opt, 'template')
    template = opt.template;
else
    template = mY; % use the mean as template
%         template = Y(:,:,1); % use the first frme as template
end

if opt.mask % create a mask, only use the image inside mask for registration
    figure, imshow(template,[])
    h = imrect;
    opt.mask_pos = h.getPosition;
    Xmin = round(opt.mask_pos(1)); 
    Ymin = round(opt.mask_pos(2));
    Xmax = Xmin + round(opt.mask_pos(3));
    Ymax = Ymin + round(opt.mask_pos(4));
    template = template(Ymin:Ymax, Xmin:Xmax);
%     h = images.roi.Rectangle(gca,'Position',[floor(para.width/2)-floor(para.height/2), 1,...
%                                              floor(para.height), floor(para.height)]); % xmin, ymin, width, height
%     title('Press Enter after the position is adjusted')
%     pause
%     Xmin = round(h.Position(1)); 
%     Ymin = round(h.Position(2));
%     Xmax = Xmin + round(h.Position(3));
%     Ymax = Ymin + round(h.Position(4));

    % only register selected 
    Zind = opt.ind_reg;
else
    Xmin = 1; Ymin = 1;
    Xmax = para.width;
    Ymax = para.height;
%     Zmin = 1;
%     Zmax = size(Y,3);
    Zind = opt.ind_reg;
end

Y = single(Y);                 % convert to single precision 
Y = Y - min(Y(:));

%% set parameters for rigid motion correction
if opt.mask
    Y_temp = Y(Ymin:Ymax, Xmin:Xmax, Zind);
else
    Y_temp = Y(:,:,Zind);
end

if strcmp(opt.mode, 'rigid')
    options_rigid = NoRMCorreSetParms(...
        'd1',           size(Y_temp,1),...
        'd2',           size(Y_temp,2),...
        'bin_width',    10,...
        'max_shift',    10,...
        'us_fac',       20,...
        'init_batch',   200,...
        'correct_bidir',0); % do not perform bi-directional scanning correction (avoid zigzag artifact)

    % perform motion correction
    tic; [Yreg_part,shifts1,template1,options_rigid] = normcorre(Y_temp, options_rigid, template); tRegistration = toc
    Yreg = Y;
    if opt.mask
        Yreg(:,:,Zind) = apply_shifts(Y(:,:,Zind), shifts1, options_rigid);
    else
        Yreg(:,:,Zind) = Yreg_part;
    end

    % save file as tiff
    % M1 = uint16(M1);
    % saveastiff(M1,'D:\=data=\80Z_imaging\img_2p\CaImAn_Results\20180503T121108_nobidir.tif');
    % saveastiff(template1,'D:\=data=\80Z_imaging\img_2p\CaImAn_Results\20180503T121108_template.tif');
end
%% set parameters for non-rigid registration
if strcmp(opt.mode, 'nonrigid')
    % clear large data in memory, to avoid out of memory problems
    options_nonrigid = NoRMCorreSetParms(...
        'd1',size(Y,1),...
        'd2',size(Y,2),...
        'grid_size',[10, 10],...
        'mot_uf',4,...
        'bin_width',200,...
        'max_shift',15,...
        'max_dev',3,...
        'us_fac',50,...
        'init_batch',200,...
        'correct_bidir',0); % do not perform bi-directional scanning correction (avoid zigzag artifact)

    % perform motion correction
    tic; [M2,shifts2,template2,options_nonrigid] = normcorre_batch(Y(Ymin:Ymax, Xmin:Xmax, Zind), options_nonrigid, template); toc
    Yreg = Y;
    if opt.mask
        Yreg(:,:,Zind) = apply_shifts(Y(:,:,Zind), shifts2, options_nonrigid);
    else
        Yreg(:,:,Zind) = M2;
    end
    % M2 = uint16(M2);
    % saveastiff(M2,'D:\=data=\80Z_imaging\img_2p\CaImAn_Results\20180503T121108_nonrigid_nonbidir.tif');
    % saveastiff(template2,'D:\=data=\80Z_imaging\img_2p\CaImAn_Results\20180503T121108_template2.tif');
end
%% compute metrics

nnY = quantile(Y(:),0.005);
mmY = quantile(Y(:),0.995);
 
[cY, mY, vY] = motion_metrics(Y,opt.bnd);
[cYreg, mYreg, vM1] = motion_metrics(Yreg,opt.bnd);
% [cM2,mM2,vM2] = motion_metrics(M2,10);
T = length(cY);
%% plot metrics
figure;
    ax1 = subplot(2,2,1); imagesc(mY,[nnY,mmY]);  axis equal; axis tight; axis off; title('mean raw data','fontsize',14,'fontweight','bold')
    ax2 = subplot(2,2,2); imagesc(mYreg,[nnY,mmY]);  axis equal; axis tight; axis off; title('mean rigid corrected','fontsize',14,'fontweight','bold')
%     ax3 = subplot(2,3,3); imagesc(mM2,[nnY,mmY]); axis equal; axis tight; axis off; title('mean non-rigid corrected','fontsize',14,'fontweight','bold')
    subplot(2,2,3); plot(1:T,cY,1:T,cYreg); legend('raw data','rigid'); title('correlation coefficients','fontsize',14,'fontweight','bold')
    subplot(2,2,4); scatter(cY,cYreg); hold on; plot([0.9*min(cY),1.05*max(cYreg)],[0.9*min(cY),1.05*max(cYreg)],'--r'); axis square;
        xlabel('raw data','fontsize',14,'fontweight','bold'); ylabel('rigid corrected','fontsize',14,'fontweight','bold');
%     subplot(2,3,6); scatter(cYreg,cM2); hold on; plot([0.9*min(cY),1.05*max(cYreg)],[0.9*min(cY),1.05*max(cYreg)],'--r'); axis square;
%         xlabel('rigid corrected','fontsize',14,'fontweight','bold'); ylabel('non-rigid corrected','fontsize',14,'fontweight','bold');
    linkaxes([ax1,ax2],'xy')
%% plot shifts        
if strcmp(opt.mode, 'rigid')
    shifts_r = squeeze(cat(3,shifts1(:).shifts));
elseif strcmp(opt.mode, 'nonrigid')
    shifts_r = squeeze(cat(3,shifts2(:).shifts));
end
% shifts_nr = cat(ndims(shifts2(1).shifts)+1,shifts2(:).shifts);
% shifts_nr = reshape(shifts_nr,[],ndims(Y)-1,T);
% shifts_x = squeeze(shifts_nr(:,1,:))';
% shifts_y = squeeze(shifts_nr(:,2,:))';

% patch_id = 1:size(shifts_x,2);
% str = strtrim(cellstr(int2str(patch_id.')));
% str = cellfun(@(x) ['patch # ',x],str,'un',0);

figure;
    ax1 = subplot(311); plot(1:T,cY,1:T,cYreg); legend('raw data','rigid'); title('correlation coefficients','fontsize',14,'fontweight','bold')
            set(gca,'Xtick',[])
    ax2 = subplot(312); plot(shifts_r(:,1),'k','linewidth',2); title('displacements along x','fontsize',14,'fontweight','bold')
            set(gca,'Xtick',[])
    ax3 = subplot(313); plot(shifts_r(:,2),'k','linewidth',2); title('displacements along y','fontsize',14,'fontweight','bold')
            xlabel('timestep','fontsize',14,'fontweight','bold')
    linkaxes([ax1,ax2,ax3],'x')

%% plot a movie with the results
% figure;
% for t = 1:1000
%     subplot(121);imagesc(Y(:,:,t),[nnY,mmY]); xlabel('raw data','fontsize',14,'fontweight','bold'); axis equal; axis tight;
%     title(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14); colormap('bone')
%     subplot(122);imagesc(Yreg(:,:,t),[nnY,mmY]); xlabel('non-rigid corrected','fontsize',14,'fontweight','bold'); axis equal; axis tight;
%     title(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14); colormap('bone')
%     set(gca,'XTick',[],'YTick',[]);
%     drawnow;
%     pause(0.01);
% end

%%
if opt.regreps % register mean images of repetitions, skip DataMat step
else
    DataMat_temp = reshape(Yreg, para.height, para.width, para.nFrame, para.nStim, length(opt.reps));
    %     DataMat_reg = reshape(Yreg, para.height, para.width, para.nFrame, para.nStim, length(opt.reps));
    [~,ind] = sort(para.order);
    tic, DataMat_temp = DataMat_temp(:,:,:,ind,:); tReorder = toc % re-arrange according to the experiment order
    tic, DataMat_temp = permute(DataMat_temp,[5, 4, 1, 2, 3]); tPermute = toc % DataMat_reg = [height, width, frames, trial, rep]
    output.DataMat_temp  = DataMat_temp;
end

output.Yreg = Yreg;
output.cYreg = cYreg;
output.mYreg = mYreg;
if strcmp(opt.mode, 'rigid')
    output.shifts = shifts1;
    output.options = options_rigid;
elseif strcmp(opt.mode, 'nonrigid')
    output.shifts = shifts2;
    output.options = options_nonrigid;
end
end