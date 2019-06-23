function ViewData(para,DataMat,opt)
% DataMat = [rep, trial, height, width, frams]

if isempty(opt.reps)
    opt.reps = 1:para.nRep;
end

if isempty(opt.tWindow)
    opt.tWindow = [para.preStim, para.preStim + para.durStim];
end

%% generate matrix of movies, for different display mode
if strcmp(opt.mode,'avgrep')
% average repetition: selected trials, averaged reps, images + videos
    nPanels = length(opt.trials);
    img_rel = zeros(nPanels,para.height,para.width);
    mov_rel = zeros(nPanels,para.height,para.width,para.nFrame);
    
    for i = 1:nPanels
        iTrial = opt.trials(i);
        mov = DataMat(opt.reps,iTrial,:,:,:);
        % average across reps
        mov_mean = squeeze(mean(mov,1));
        % calculate deltaF/F (movie)
        img_base = squeeze(mean(mov_mean(:,:,1:floor(para.fr*para.preStim)),3)); % pre-stimulus: baseline        
        img_base = repmat(img_base,1,1,para.nFrame);
        mov_rel(i,:,:,:) = (mov_mean - img_base)./img_base;
        % calculate deltaF/F (averaged image)
        img_rel(i,:,:) = squeeze(mean(mov_rel(i,:,:,floor(para.fr*opt.tWindow(1))+1 : floor(para.fr*opt.tWindow(2))),4));  
    end
    
    fnametemp = para.filename;
    vnametemp = [para.pathname, para.filename(1:end-4), '_trial', ...
        num2str(opt.trials(1)),'-',num2str(opt.trials(end)),'.avi']; 


elseif strcmp(opt.mode,'allrep')
    nPanels = length(opt.reps) + 1;
    if length(opt.trials) == 1
        iTrial = opt.trials; 
    else
        disp('Error: number of trial selected must be 1')
    end

    % ==== calculate relative amplitude change ====
    mov = squeeze(DataMat(opt.reps,iTrial,:,:,:)); 
    mov_rel = mov; % #Rep x height x width x frames
    img_rel = zeros(length(opt.reps)+1, para.height, para.width); % response mean image for each rep; the last one is averaged
    for i = 1:length(opt.reps)
        mov_temp = squeeze(mov(i,:,:,:)); % 75 x 120 x 25
        img_base = squeeze(mean(mov_temp(:,:,1:floor(para.fr*para.preStim)),3));
        img_base = repmat(img_base,1,1,para.nFrame);
        mov_rel(i,:,:,:) = (mov_temp - img_base)./img_base;
        img_rel(i,:,:) = mean(mov_rel(i,:,:,floor(para.fr*opt.tWindow(1))+1 : floor(para.fr*opt.tWindow(2))),4);
    end
    
     % ==== calculate mean movie, put it after all reps ==== 
    mov_mean = squeeze(mean(mov)); %averaged across reps, height x width x frames
    img_base = squeeze(mean(mov_mean(:,:,1:floor(para.fr*para.preStim)),3));
    img_base = repmat(img_base,1,1,para.nFrame);
    mov_rel(i+1,:,:,:) = (mov_mean - img_base)./img_base;
    img_rel(i+1,:,:) = mean(mov_rel(i+1,:,:,floor(para.fr*opt.tWindow(1))+1 : floor(para.fr*opt.tWindow(2))),4);
    
    fnametemp = para.filename;
    vnametemp = [para.pathname, para.filename(1:end-4), '_trial', num2str(iTrial),'.avi'];
end


%% save video, start a video writer object
switch opt.plotMode
    case 'combined'
    if opt.saveON
        if opt.soundON % video with sound
            [file,path] = uigetfile('*.wav','Select the sound file','\\FANTASIA-DS3617\Test_Imaging\=Sounds=');
            [sounddata,fs] = audioread(fullfile(path,file)); 
            videoFWriter = vision.VideoFileWriter(vnametemp,...
            'FileFormat',       'AVI',...
            'AudioInputPort',   true,...
            'FrameRate',        para.fr,...
            'VideoCompressor',	'None (uncompressed)',...
            'AudioDataType',    'int16');
            SoundBatchSampleNum =   round(fs/para.fr);
            SoundSeq =              sounddata;
        else % video without sound
            videoFWriter = vision.VideoFileWriter(vnametemp,...
            'FileFormat',       'AVI',...
            'AudioInputPort',   false,...
            'FrameRate',        para.fr,...
            'VideoCompressor',	'None (uncompressed)',...
            'AudioDataType',    'int16');
        end
    end

    % arrange subplots 
    [p,n] = numSubplots(nPanels);
    % construct a combined matrix for mean images
    img_rel = permute(img_rel,[2,3,1]); % height x width x trials
    if nPanels < n % pad with 0 if one block is left empty
        img_rel(:,:,nPanels+1:n) = zeros(para.height,para.width,n-nPanels);
    end
    % rearrange the matrix above into a large matrix
    img_all = [];    
    for k = 1:p(1)
        img_all = [img_all;reshape(img_rel(:,:,p(2)*(k-1)+1:p(2)*k),para.height,para.width*p(2))];
    end

    figure,
    % display the averaged reponse first
    h = imagesc(img_all,opt.ampLimit); colormap('jet'); colorbar; axis image
    pause
    for i = 1:para.nFrame
        mov_temp = mov_rel(:,:,:,i);
        % construct a combined matrix for each frame
        mov_temp = permute(mov_temp,[2,3,1]);
        if nPanels < n % pad with 0 if one block is left empty
            mov_temp(:,:,nPanels+1:n) = zeros(para.height,para.width,n-nPanels);
        end
        mov_all = [];
        for k = 1:p(1)
            mov_all = [mov_all;reshape(mov_temp(:,:,p(2)*(k-1)+1:p(2)*k),para.height,para.width*p(2))];
        end

        set(h,'CData',mov_all)
        title(['time = ',num2str(i*0.2,'%-5.1f')])
        pause(1/para.fr)
    %         pause
        if opt.saveON    
            % save video and audio
            MAX =           opt.ampLimit(2);
            MIN =           min(min(mov_all));
            mov_all =       (2^8-1).* (mov_all - MIN)./(MAX - MIN);
            mov_all =       uint8(mov_all);
            frame =         repmat(mov_all,1,1,3);
            frame(:,:,1) =  mov_all*0; % use green color
            frame(:,:,3) =  mov_all*0;
            if opt.soundON
                if i*SoundBatchSampleNum > length(SoundSeq) % if end of current frame is longer than sound (by <1 segment)
                    videoFWriter(frame,SoundSeq((i-1)*SoundBatchSampleNum+1:end) );
                elseif (i-1)*SoundBatchSampleNum+1 > length(SoundSeq) % if 
                    videoFWriter(frame);
                else            
                    videoFWriter(frame,SoundSeq((i-1)*SoundBatchSampleNum+1:i*SoundBatchSampleNum) );
                end
            else
                videoFWriter(frame);
            end
        end
    end
    set(h,'CData',img_all)
    if opt.saveON % write the last frame and release object videoFWriter
        MAX =           opt.ampLimit(2);
        MIN =           min(min(img_all));
        img_all =       (2^8-1).*(img_all - MIN)./(MAX - MIN);
        img_all =       uint8(img_all);
        frame =         repmat(img_all,1,1,3);
        frame(:,:,1) =  img_all*0;
        frame(:,:,3) =  img_all*0;
        videoFWriter(frame);
        release(videoFWriter) 
    end
 %% plot separately, cannot save videos
    case 'separate'
    % show movie as separate matrices (the last is the averaged across reps)
    switch opt.mode
        case 'allrep'
        [p,n] = numSubplots(nPanels);
        figure,
        for iRep = 1:nPanels
            temp = squeeze(mov_rel(iRep,:,:,:));
            subplot(p(1),p(2),iRep); 
            h(iRep) = imagesc(temp(:,:,1),opt.ampLimit);colorbar;
            axis image
            if iRep == nPanels
                title('Averaged across reps')
            else
                title(['Rep.Num.',num2str(opt.reps(iRep))])
            end
        end
        pause
        for i = 1:para.nFrame
            for iRep = 1:nPanels
            set(h(iRep),'CData',squeeze(mov_rel(iRep,:,:,i)))
            end
            pause(1/para.fr)
        end
        
        case 'avgrep'
        [p,n] = numSubplots(nPanels);
        figure,
        for iTrial = 1:nPanels
            temp = squeeze(mov_rel(iTrial,:,:,:));
            subplot(p(1),p(2),iTrial); 
            h(iTrial) = imagesc(temp(:,:,1),opt.ampLimit);colorbar;
            axis image
            title(['Rep.Num.',num2str( opt.trials(iTrial) )])
        end
        pause
        for i = 1:para.nFrame
            for iTrial = 1:nPanels
            set(h(iTrial),'CData',squeeze(mov_rel(iTrial,:,:,i)))
            end
            pause(1/para.fr)
        end    
        
    end
end
    
end


% ===========================================
% reference code is from Xindong Song (for saving video with audio)
% ===========================================
% for saving through VideoFWriter System object (2015/08/17)
% fnametemp = F.FileName{1};
% vnametemp = [F.PathName, F.FileName{1}(1:end-4), '.avi'];
% sounddata = S.SesSoundWave; 
% framerate = 80/P.ProcBinFrame;
% T.hVideoFileWriter = vision.VideoFileWriter(vnametemp,...
%     'FileFormat',       'AVI',...
%     'AudioInputPort',   true,...
%     'FrameRate',        framerate,...
%     'VideoCompressor',	'None (uncompressed)',...
%     'AudioDataType',    'int16');
% 
% SoundBatchSampleNum =   round(100e3/framerate);
% SoundSeq =              sounddata;
% 
% for i = 1:R.N_Fpc
%     disp(num2str(i));
% 	T.Frame =               squeeze(T.Polarity *R.FpcPhPw_NormCycMean(i,:,:)/T.saturation);
% 	T.Frame =               max(max(T.Frame, 0), 0);
% 	T.Frame =               min(min(T.Frame, 1), 1);  
%     T.FrameInt =            uint8(255*T.Frame);
%     T.FrameOut(:,:,1) =     T.FrameInt*0;
%     T.FrameOut(:,:,2) =     T.FrameInt;
%     T.FrameOut(:,:,3) =     T.FrameInt*0;
%     T.FrameOut =            uint8(T.FrameOut);
%     
%     step(   T.hVideoFileWriter,...
%             T.FrameOut,...
%             SoundSeq((i-1)*SoundBatchSampleNum+1:i*SoundBatchSampleNum)  );
%     imagesc(    T.Frame);
%     title(      sprintf('%5.1f s', i/(80/Xin.P.ProcBinFrame)) );   
%     pause(0.1);
% end
% 
% release(T.hVideoFileWriter);