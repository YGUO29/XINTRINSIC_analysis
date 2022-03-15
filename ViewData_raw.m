function X = ViewData_raw(DataMat, para, opt)
% for viewing a continous stimulus session, baseline set as the average
% amplitude of the entire session

% DataMat = [height, width, frams]
mov = squeeze(DataMat); % eliminate trial and/or repetition dimensions
if length(size(mov)) == 4 && strcmp(opt.mode, 'avgrep') % more than one reps, average across reps
    mov = squeeze(mean(mov,1));
elseif length(size(mov)) == 4 && strcmp(opt.mode, 'allrep')
    mov = permute(mov, [2 3 4 1]);
    mov = reshape(mov, para.height, para.width, size(mov, 3)*size(mov, 4));
    para.nFrame = size(mov,3); 
else % only one rep
end
    
baseline    = repmat(mean(mov,3), 1, 1, para.nFrame);
deltaF      = mov - baseline;
% tic, 
mov_rel     = deltaF./baseline;
X           = reshape(mov_rel, para.height*para.width, size(mov_rel,3));
X = X';
% toc

% for i = 1:para.nFrame
%     imagesc(mov_rel(:,:,i)); colormap('gray'); colorbar; axis image
%     pause
% end
% Delta   = bsxfun(@minus, mov, mean(mov,3));
% mov_rel = bsxfun(@rdivide, Delta, mean(mov,3));

%% save video, start a video writer object
if opt.saveON
    fnametemp = para.filename;
    vnametemp = [cd,'\', para.filename(1:end-4), '.avi'];
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

   
    figurex,
    % display the averaged reponse first
    h = imagesc(mean(mov_rel,3),opt.ampLimit); 
    colormap('jet'); colorbar; 
%     CT = cbrewer('div', 'RdBu', 255);
%     colormap(CT); colorbar
    axis image
    pause
    for i = 1:para.nFrame
        mov_temp = mov_rel(:,:,i);
        set(h,'CData',mov_temp)
        title(['time = ',num2str(i*0.2,'%-5.1f')])
        pause(1/para.fr/5)
        if opt.saveON    
            % save video and audio
            MAX =           opt.ampLimit(2);
%             MIN =           min(min(mov_temp));
            MIN =           opt.ampLimit(1);
            mov_temp =       (2^8-1).* (mov_temp - MIN)./(MAX - MIN);
            mov_temp =       uint8(mov_temp);
            frame =         repmat(mov_temp,1,1,3);
         
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
    if opt.saveON % write the last frame and release object videoFWriter
        release(videoFWriter) 
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