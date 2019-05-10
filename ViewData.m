function ViewData(para,DataMat,opt)
% DataMat = [rep, trial, height, width, frams]

% mode1: selected trials, averaged reps, images + videos
% mode2: 1 trial, all reps separate, images + movies
% mode3: 1 trial, averaged reps, movie
%%
if strcmp(opt.mode,'avgrep')
% mode1: selected trials, averaged reps, images + videos

    nTrials = length(opt.trials);
    [p,n] = numSubplots(nTrials);
    figure,
    for i = 1:nTrials
        iTrial = trials(i);
        mov = DataMat(:,iTrial,:,:,:);
        % average across reps
        mov_mean = squeeze(mean(mov,1));
        % calculate deltaF/F (movie)
        img_base = squeeze(mean(mov_mean(:,:,1:floor(para.fr*para.preStim)),3)); % pre-stimulus: baseline        
        img_base = repmat(img_base,1,1,para.nFrame);
        mov_rel = (mov_mean - img_base)./img_base;
        % calculate deltaF/F (averaged image)
        img_rel = squeeze(mean(mov_rel(:,:,floor(para.fr*para.preStim)+1 : floor(para.fr*(para.preStim + para.durStim))),3));

        subplot(p(1),p(2),i)
        imagesc(img_rel,[0,opt.ampLimit])
        axis image
        colorbar
    end

%%    
elseif strcmp(opt.mode,'allrep')
% mode2: 1 trial, all reps separate, images + movies
    if length(opt.trials) == 1
        iTrial = opt.trials; 
    else
        disp('Error: number of trial selected must be 1')
    end

    % ==== calculate relative amplitude change ====
    mov = squeeze(DataMat(:,iTrial,:,:,:)); 
    mov_rel = mov;
    img_rel = zeros(para.nRep+1, para.height, para.width); % for each rep, the last one is averaged
    for i = 1:para.nRep
        mov_temp = squeeze(mov(i,:,:,:)); % 75 x 120 x 25
        img_base = squeeze(mean(mov_temp(:,:,1:floor(para.fr*para.preStim)),3));
        img_base = repmat(img_base,1,1,para.nFrame);
        mov_rel(i,:,:,:) = (mov_temp - img_base)./img_base;
        img_rel(i,:,:) = mean(mov_rel(i,:,:,floor(para.fr*para.preStim)+1 : floor(para.fr*(para.preStim + para.durStim))),4);
    end

    % ==== calculate mean movie ==== 
    mov_mean = squeeze(mean(mov));
    img_base = squeeze(mean(mov_mean(:,:,1:floor(para.fr*para.preStim)),3));
    img_base = repmat(img_base,1,1,para.nFrame);
    mov_rel(i+1,:,:,:) = (mov_mean - img_base)./img_base;
    img_rel(i+1,:,:) = mean(mov_rel(i+1,:,:,floor(para.fr*para.preStim)+1 : floor(para.fr*(para.preStim + para.durStim))),4);

    if strcmp(opt.plotMode,'combined')
    % ==== show movie as combined matrix (the last is the averaged across reps)    
    
    % for saving through VideoFWriter System object (2015/08/17)
    if opt.saveON    
        fnametemp = para.filename;
        vnametemp = [para.pathname, para.filename(1:end-4), '.avi'];
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
    end

        % convert separate matrices to a single matrix    
        [p,n] = numSubplots(para.nRep+1); % decide how to arrange them
        % construct a combined matrix for mean images
        img_rel = permute(img_rel,[2,3,1]);
        if para.nRep+1 < n % pad with 0 if one block is left empty
            img_rel(:,:,para.nRep+2:n) = zeros(para.height,para.width,n-para.nRep-1);
        end
        img_all = [];    
        for k = 1:p(1)
            img_all = [img_all;reshape(img_rel(:,:,p(2)*(k-1)+1:p(2)*k),para.height,para.width*p(2))];
        end

        figure,
        h = imagesc(img_all,[0,opt.ampLimit]); colorbar; axis image
        pause
        for i = 1:para.nFrame
            mov_temp = mov_rel(:,:,:,i);
            % construct a combined matrix for each frame
            mov_temp = permute(mov_temp,[2,3,1]);
            if para.nRep+1 < n % pad with 0 if one block is left empty
            mov_temp(:,:,para.nRep+2:n) = zeros(para.height,para.width,n-para.nRep-1);
            end
            mov_all = [];
            for k = 1:p(1)
                mov_all = [mov_all;reshape(mov_temp(:,:,p(2)*(k-1)+1:p(2)*k),para.height,para.width*p(2))];
            end

            set(h,'CData',mov_all)
            title(['time = ',num2str(i*0.2,'%-5.1f')])
            pause(0.2)
            
            % save video and audio
            MAX =           opt.ampLimit;
            MIN =           min(min(mov_all));
            mov_all =       (2^8-1).*(mov_all - MIN)./(MAX - MIN);
            frame =         repmat(mov_all,3,1);
            frame(:,:,1) =  mov_all*0;
            frame(:,:,3) =  mov_all*0;
            Filename =      videoFWriter(frame,SoundSeq((i-1)*SoundBatchSampleNum+1:i*SoundBatchSampleNum) );
        end
        set(h,'CData',img_all)
        release(videoFWriter)
    end

    if strcmp(opt.plotMode,'separate')
    % ==== show movie as separate matrices (the last is the averaged across reps)
        [p,n] = numSubplots(para.nRep+1);
        figure,
        for iRep = 1:para.nRep+1
            temp = squeeze(mov_rel(iRep,:,:,:));
            subplot(p(1),p(2),iRep); 
            h(iRep) = imagesc(temp(:,:,1),[0,opt.ampLimit]);colorbar;
            axis image
            if iRep == paara.nRep+1
                title('Averaged across reps')
            else
                title(['Rep.Num.',num2str(iRep)])
            end
        end
        pause
        for i = 1:para.nFrame
            for iRep = 1:para.nRep+1
            set(h(iRep),'CData',squeeze(mov_rel(iRep,:,:,i)))
            end
            pause(0.2)
        end
    end
    
end




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