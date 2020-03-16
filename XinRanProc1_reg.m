%% Xintrinsic Processing 1 - with rigid registration (added by Yueqi Oct. 2019)


%% DATA BINNNING

%% Locate Files 
clear all
[F.FileName, F.PathName, F.FilterIndex] = uigetfile(...
    'Z:\2018.03 T2 (Marmoset 80Z, Xintrinsic, Green & Fluo)\*.rec',...
    'Select raw recording files to process',...
    'MultiSelect',              'On');
if F.FilterIndex == 0
    clear F;                    % nothing selected
    return
end
if iscell(F.FileName) == 0      % single file selected
    F.FileName = {F.FileName};
end

disp(['Xintrinsic Processing Stage 1 (spatiotemporal binning) is about to start on ' ...
    num2str(length(F.FileName)) ' files']);

%% Parmateters
% R: Recorded
% P: Processed
% T: Temporal
P.ProcMask = 0; % load a surface picture to draw mask? (for registration)

R.SysCamPixelHeight =           300;
R.SysCamPixelWidth =            480;
R.SysCamPixelBinNum =           4;
R.SysCamFrameRate =             80;

% Bin 16 x 16 pixels together   (F) 
P.ProcPixelBinNum =             4; %1: original resolution for registration
P.ProcPixelHeight =             R.SysCamPixelHeight/P.ProcPixelBinNum;
P.ProcPixelWidth =              R.SysCamPixelWidth/P.ProcPixelBinNum;

% Bin 16 frames together        (Proc)
P.ProcFrameRate =               5;
P.ProcFrameBinNum =             R.SysCamFrameRate/P.ProcFrameRate; 

T.hWaitbar =                    waitbar(0, 'processing');

%% DATA BINNING
for i = 1: length(F.FileName)
   
    %% Load 'S'
    T.filename = [F.PathName, F.FileName{i}];
    load([T.filename(1:end-3) 'mat']);  
    disp([  'Processing: "', F.FileName{i}, ...
            '" with the sound: "', S.SesSoundFile, '"']);
    
    %% Parameter initialization for Spatial & Temporal Binning
    
    R.SesTrlNumTotal =          length(S.SesTrlOrderVec);
    R.SysCamFramePerTrial =     S.TrlDurTotal * R.SysCamFrameRate;
    P.ProcFramePerTrial =       S.TrlDurTotal * P.ProcFrameRate;   
    P.ProcFrameNumTotal =       S.SesFrameTotal / P.ProcFrameBinNum;
    
    T.fid =                     fopen(T.filename);
    
	P.RawMeanPixel =            zeros(1, S.SesFrameTotal);
    P.RawMeanPower =            zeros(1, S.SesFrameTotal);                                
    P.ProcMeanPixel =           zeros(1, P.ProcFrameNumTotal);
    P.ProcMeanPower =           zeros(1, P.ProcFrameNumTotal);
    
    P.ProcDataMat =             zeros(...
                                    S.SesCycleNumTotal,...
                                    S.TrlNumTotal,...
                                    P.ProcPixelHeight,...
                                    P.ProcPixelWidth,...
                                    P.ProcFramePerTrial...
                                    );                         
    for j = 1:S.SesCycleNumTotal
        for k = 1:S.TrlNumTotal
            m = (j-1)*S.TrlNumTotal + k;
            %% Update GUI
            waitbar(m/R.SesTrlNumTotal, T.hWaitbar,...
                ['finishing ',...
                sprintf('%d out of %d total trials in the session',...
                    m, R.SesTrlNumTotal)] );       
            %% Read Data Batch    
            T.DataRaw = 	fread(T.fid, [...
                R.SysCamPixelHeight * R.SysCamPixelWidth, ...
                R.SysCamFramePerTrial],...
                'uint16');
            %% Frame #, Trial order # location        
            T.RecFramesCurrent =    ((m-1)*	R.SysCamFramePerTrial +1):...
                                    (m*     R.SysCamFramePerTrial);
            T.ProcFramesCurrent =   ((m-1)*	P.ProcFramePerTrial +1):...
                                    (m*     P.ProcFramePerTrial);        
            T.TrialOrder =          S.SesTrlOrderVec(m);            
            %% Image Processing
            T.PixelMeanRaw =        mean(T.DataRaw, 1);
            T.PixelMeanBinned =     mean( reshape(...
                                                T.PixelMeanRaw,...
                                                P.ProcFrameBinNum,...
                                                P.ProcFramePerTrial), 1 );
            P.RawMeanPixel(T.RecFramesCurrent) =    T.PixelMeanRaw;
            P.ProcMeanPixel(T.ProcFramesCurrent) =  T.PixelMeanBinned;  
            
            T.ImageS1 =         reshape(T.DataRaw,...  
                P.ProcPixelBinNum,     P.ProcPixelHeight, ...
                P.ProcPixelBinNum,     P.ProcPixelWidth, ...
                P.ProcFrameBinNum,     P.ProcFramePerTrial); 
            T.ImageS2 =         sum(T.ImageS1, 1);  
            T.ImageS3 =         sum(T.ImageS2, 3); 
            T.ImageS4 =         sum(T.ImageS3, 5);
            T.ImageS5 =         squeeze(T.ImageS4);
            P.ProcDataMat(j, k, :, :, :) =...
                                        T.ImageS5;  % in the order of experiment                       
        end    
    end
    
    %% Registration with NoRMCorre
    tic, P.ProcDataMat = permute(P.ProcDataMat,[3,4,5,2,1]); time.Permute1 = toc; % from [rep, trial, h, w, f] to [h, w, f, trial, rep]
    Y = reshape(P.ProcDataMat,...
        P.ProcPixelHeight,...
        P.ProcPixelWidth,...
        P.ProcFramePerTrial*S.TrlNumTotal*S.SesCycleNumTotal); 
    template = Y(:,:,1); % use the first frme as template
    if P.ProcMask
%         [F.FileNamePic, F.PathNamePic, F.FilterIndexPic] = uigetfile(...
%         [F.PathName,'*.tif'],...
%         'Select a surface image',...
%         'MultiSelect', 'off');
%         I = imread(fullfile(F.PathNamePic,F.FileNamePic));
%         I = double(imresize(I,[P.ProcPixelHeight,P.ProcPixelWidth]));
    
        figure, imshow(template,[])
        % h = images.roi.Circle(gca,'Center',[floor(para.width/2) floor(para.height/2)],'Radius',floor(para.height/2)); 
        h = images.roi.Rectangle(gca,'Position',[floor(P.ProcPixelWidth/2)-floor(P.ProcPixelHeight/2), 1,...
                                                 floor(P.ProcPixelHeight), floor(P.ProcPixelHeight)]); % xmin, ymin, width, height
        % h = images.roi.Ellipse(gca,'Center',[floor(para.width/2) floor(para.height/2)],'Semiaxes',[40 20]); 
        % h = images.roi.Polygon(gca,'Position',[1 1; 1 para.height - 30; 30 para.height; para.width para.height; para.width 1]); 
        title('Press Enter after the position is adjusted')
        pause
        Xmin = round(h.Position(1));    Xmax = Xmin + round(h.Position(3));
        Ymin = round(h.Position(2));    Ymax = Ymin + round(h.Position(4));
        Zmin = 1;                       Zmax = size(Y,3);
        
    else
        Xmin = 1; Xmax = P.ProcPixelWidth;
        Ymin = 1; Ymax = P.ProcPixelHeight;
        Zmin = 1; Zmax = size(Y,3);
    end
    
    Y = single(Y);                 % convert to single precision 
%     T = size(Y,ndims(Y));
    Y = Y - min(Y(:));
    
    % set parameters for motion correction
   
    options_rigid = NoRMCorreSetParms(...
        'd1',           size(Y,1),...
        'd2',           size(Y,2),...
        'bin_width',    10,...
        'max_shift',    20,...
        'us_fac',       20,...
        'init_batch',   200,...
        'correct_bidir',0); % do not perform bi-directional scanning correction (avoid zigzag artifact)

    % perform motion correction
    disp('Registration start...')
    tic; [Yreg_part,shifts1,template1,options_rigid] = normcorre(Y(Ymin:Ymax, Xmin:Xmax, Zmin:Zmax),options_rigid, template); time.Registration = toc
    if P.ProcMask
        Yreg = apply_shifts(Y,shifts1,options_rigid);
    else
        Yreg = Yreg_part;
    end
    disp(['Registration finished, time elapsed = ',num2str(time.Registration)])
    P.ProcDataMat = reshape(Yreg, P.ProcPixelHeight, P.ProcPixelWidth, P.ProcFramePerTrial, S.TrlNumTotal, S.SesCycleNumTotal);
    for iSes = 1:S.SesCycleNumTotal
        [~,ind] = sort(S.SesTrlOrderMat(iSes,:));
        tic, P.ProcDataMat(:,:,:,:,iSes) = P.ProcDataMat(:,:,:,ind,iSes); time.Reorder2 = toc % re-arrange according to the experiment order
    end
    tic, P.ProcDataMat = permute(P.ProcDataMat,[5, 4, 1, 2, 3]); time.Permute2 = toc % DataMat_reg = [height, width, frames, trial, rep]

    
        %% Power Processing
%         P.RawMeanPower =    mean(S.SesPowerMeter, 2)';
%         P.ProcMeanPower =   mean(reshape(P.RawMeanPower,...
%                                     P.ProcFrameBinNum,...
%                                     P.ProcFrameNumTotal), 1 );
    %% Show Figure
    T.timeraw =                 (1:S.SesFrameTotal)/R.SysCamFrameRate;
    T.timebinned =              (1:P.ProcFrameNumTotal)/P.ProcFrameRate;
    figure(     'Name',         F.FileName{i});
    subplot(2,1,1);
    [T.hAx,~,~] =       plotyy(	T.timeraw,      P.RawMeanPower, ...
                                T.timeraw,      P.RawMeanPixel);
    xlabel(T.hAx(1),        'Time (sec)');
    ylabel(T.hAx(1),        'Power Mean (volt)');
    ylabel(T.hAx(2),        'Pixel Mean (ADU)');
    subplot(2,1,2);
    [T.hAx, T.hP1, T.hP2] = ...
                        plotyy(	T.timebinned,   P.ProcMeanPower, ...
                                T.timebinned,   P.ProcMeanPixel);    
%     T.hP2.LineWidth =       2;                     
    xlabel(T.hAx(1),        'Time (sec)');
    ylabel(T.hAx(1),        'Power Mean (volt)');
    ylabel(T.hAx(2),        'Pixel Mean (ADU)');
    if P.ProcPixelBinNum == 1
        save([T.filename(1:end-4) '_P1_reg_300x480.mat'], 'P', '-v7.3');     
    elseif P.ProcPixelBinNum == 4
        save([T.filename(1:end-4) '_P1_reg.mat'], 'P', '-v7.3'); 
    end
end

close(T.hWaitbar);
disp('All files are processed');
return;
