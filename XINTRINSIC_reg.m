% registration with matlab function
function XINTRINSIC_reg(DataMat, I, para, opt)

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
%%
Yreg = Y;
[cY,mY,vY] = motion_metrics(Y,5);

fixed = Y(:,:,1);

ind = find(cY < 0.9999);
nRegFrames = length(ind);
moving = Y(:,:,ind);
% movingRegistered = moving;
[optimizer,metric] = imregconfig('monomodal');

tic
for i = 1:nRegFrames
    tform = imregtform(moving(:,:,i),fixed,'rigid', optimizer, metric);
    Yreg(:,:,ind(i)) = imwarp(moving(:,:,i),tform,'OutputView',imref2d(size(fixed)));
    i
end
toc



[cM1,mM1,vM1] = motion_metrics(Yreg,5);

DataMat_reg = reshape(Yreg, para.height, para.width, para.nFrame, para.nStim, para.nRep);
[~,ind] = sort(para.order);
tic, DataMat_reg = DataMat_reg(:,:,:,ind,:); tReorder = toc % re-arrange according to the experiment order
tic, DataMat_reg = permute(DataMat_reg,[5, 4, 1, 2, 3]); tPermute = toc % DataMat_reg = [height, width, frames, trial, rep]

end