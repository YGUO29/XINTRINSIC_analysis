%% Xintrinsic Processing 1 
%% DATA BINNNING

%% Locate Files 
clear all
[F.FileName, F.PathName, F.FilterIndex] = uigetfile(...
    '\\FANTASIA-DS3617\Test_Imaging\2018.03 T2 (Marmoset 80Z, Xintrinsic, Green & Fluo)\*.rec',...
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

R.SysCamPixelHeight =           300;
R.SysCamPixelWidth =            480;
R.SysCamPixelBinNum =           4;
R.SysCamFrameRate =             80;

% Bin 16 x 16 pixels together   (F)
P.ProcPixelBinNum =             1;
P.ProcPixelHeight =             R.SysCamPixelHeight/P.ProcPixelBinNum;
P.ProcPixelWidth =              R.SysCamPixelWidth/P.ProcPixelBinNum;

% Bin 16 frames together        (Proc)
P.ProcFrameRate =               20;
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
    
    R.SesTrlNumTotal =          S.TrlNumberTotal;
    R.SysCamFramePerTrial =     S.TrlDurTotal * R.SysCamFrameRate;
    P.ProcFramePerTrial =       S.TrlDurTotal * P.ProcFrameRate;   
    P.ProcFrameNumTotal =       S.SesFrameTotal / P.ProcFrameBinNum;
    
    T.fid =                     fopen(T.filename);
    
	P.RawMeanPixel =            zeros(1, S.SesFrameTotal);
    P.RawMeanPower =            zeros(1, S.SesFrameTotal);                                
    P.ProcMeanPixel =           zeros(1, P.ProcFrameNumTotal);
    P.ProcMeanPower =           zeros(1, P.ProcFrameNumTotal);
    
    P.ProcDataMat =             zeros(...
                                    S.SesCycleTotal,...
                                    S.TrlNumberTotal,...
                                    P.ProcPixelHeight,...
                                    P.ProcPixelWidth,...
                                    P.ProcFramePerTrial...
                                    );                         
    for j = 1:S.SesCycleTotal
        for k = 1:S.TrlNumberTotal
            m = (j-1)*S.TrlNumberTotal + k;
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
%             T.TrialOrder =          S.SesTrlOrderVec(m); 
            T.TrialOrder =          1;            

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
            P.ProcDataMat(j, T.TrialOrder, :, :, :) =...
                                        T.ImageS5;                         
        end    
    end
        %% Power Processing
        P.RawMeanPower =    mean(S.SesPowerMeter, 2)';
        P.ProcMeanPower =   mean(reshape(P.RawMeanPower,...
                                    P.ProcFrameBinNum,...
                                    P.ProcFrameNumTotal), 1 );
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
%     save([T.filename(1:end-4) '_P1.mat'], 'P', '-v7.3');     
    T.SesName = strsplit(T.filename,'\');
    T.SesName = T.SesName{end};
    T.SesName = T.SesName(1:end-4);
    save(['\\FANTASIA-DS3617\Proc_Imaging_Yueqi\Marmoset_imaging\data_wf\80Z\',T.SesName, '_P1_150x240@20fps.mat'], 'P', '-v7.3');     
end

close(T.hWaitbar);
disp('All files are processed');
return;
