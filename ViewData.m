function [X, mov_rel, mov_rel_sep] = ViewData(DataMat, para, opt)
% convert pre-processed data to data matrix for later process
% DataMat = [rep, trial, height, width, frams]

% demo usage:
%     opt.ampLimit      = 0.003 .*[-1 1];
%     opt.reps          = [1:11+3, 11+6:para.nRep];
%     opt.trials        = 6; 
%     opt.omit_trials   = [12 16 30 36 48 48 16 16];
%     opt.omit_reps     = [1 1 1 1 1 5 7 10]; 
%     opt.tWindow       = [para.preStim + 4, para.preStim + para.durStim]; % intrinsic natural sound: integrate until 4s after sound offset
%     opt.p             = [8, 6]; % subplot rows and columes; 
%     opt.mode          = 'allrep'; % allrep or avgrep
%     opt.plotMode      = 'combined'; % combined or separate
%     opt.saveON        = 1;
%     opt.soundON       = 1;
%     para.pathname     = 'D:\SynologyDrive\=data=\Marmoset_imaging\video_wf\';
%     figurex;
%     [X, DataMat_norm] = ViewData(DataMat, para, opt); % X may contain NaNs if there are masked pixels
    
% ======================= initialization ===================================
if ~isfield(opt, 'reps')
    opt.reps = 1:para.nRep; % default: all repetitions
end

% ==== window for averaging response ====
if ~isfield(opt, 'tWindow') % not specified, make it equal to stimulus duration
    opt.tWindow = [para.preStim, para.preStim + para.durStim];
    opt.tWindow = repmat(opt.tWindow, para.nStim, 1);
elseif isfield(opt, 'tWindow') && length(opt.tWindow) == 2 % specified, same durations
    opt.tWindow = repmat(opt.tWindow, para.nStim, 1);    
else % specified, different durations
end

if ~isfield(opt, 'trials')
    opt.trials = 1:para.nStim; % default: all trials
end
if ~isfield(opt, 'ampLimit')
    opt.ampLimit = 0.1.*[-1 1]; % default: +/- 10%
end
if ~isfield(opt, 'mode')
    opt.mode = 'avgrep'; % 'avgrep' (averaged across repetitions) or 'allrep' (all repetitions plotted separately)
end
if ~isfield(opt, 'plotMode')
    opt.plotMode    =  'combined'; % plot figures as 'combined' or 'separate' (video saving is only available for 'combined' mode)
end
if ~isfield(opt, 'saveON')
    opt.saveON      = 0; % save as video
end
if ~isfield(opt, 'soundON')
    opt.soundON     = 0; % save the sound for this recording as well
end

if ~isfield(opt, 'color') % default - save and display in jet
    opt.color     = 'jet'; % save the sound for this recording as well
    cmap          = jet(256);
else
    if ischar(opt.color)
        switch opt.color 
            case 'green'; cmap = jet(256); % save in green mode, display in jet
            otherwise; cmap = eval([opt.color, '(256)']); cmap = flipud(cmap); % save and display in specified cmap 
        end
    else % save and display in cbrewer maps
        if strcmp(opt.color{2}, 'RdBu'); cmap = flipud(cbrewer(opt.color{1}, opt.color{2}, 256)); % flip red and blue for RdBu map (red = positive)
        else cmap = cbrewer(opt.color{1}, opt.color{2}, 256); end
    end
end


%% throw out moving trials
% lengths of these two fields need to be matched
if isfield(opt, 'omit_trials') || isfield(opt, 'omit_reps')
    for i = 1:length(opt.omit_trials)
        DataMat(opt.omit_reps(i), opt.omit_trials(i), :, :, :) = NaN;
    end
end


%% Generate matrix of movies, for different display mode
%% Averaged across repetitions (across opt.reps)
if strcmp(opt.mode,'avgrep')
    nPanels = length(opt.trials);
    img_rel = zeros(nPanels,para.height,para.width); % img = mean response, rel = relative values
    mov_rel = zeros(nPanels,para.height,para.width,para.nFrame);  % mov = frames of the trial, rel = relative values
    mov_rel_sep = zeros(length(opt.reps), nPanels,para.height,para.width,para.nFrame);
    % average all reps
%     mov = squeeze(mean(DataMat(opt.reps,:,:,:,:),1)); % trials x height x width x frames

    for i = 1:nPanels
        iTrial = opt.trials(i);
%         % get trial iTrial, all repetitions
%         mov_mean = squeeze(mov(iTrial, :, :, :));            
%         % calculate deltaF/F (movie)
%         img_base = squeeze(mean(mov_mean(:,:,1:floor(para.fr*para.preStim)), 3, 'omitnan')); % pre-stimulus: baseline
%         img_base = repmat(img_base,1,1,para.nFrame);
%         mov_rel(i,:,:,:) = (mov_mean - img_base)./img_base;
%         % calculate deltaF/F (averaged image, time window determined by opt.tWindow)
%         img_rel(i,:,:) = squeeze(mean(mov_rel(i,:,:,floor(para.fr*opt.tWindow(iTrial, 1))+1 : floor(para.fr*opt.tWindow(iTrial, 2))),4, 'omitnan'));  
        for j = 1:length(opt.reps)
            iRep = opt.reps(j);
            img_base =  squeeze(mean(DataMat(iRep, iTrial, :, :, 1:floor(para.fr*para.preStim)), 5, 'omitnan')); 
            img_base = repmat(img_base, 1, 1, para.nFrame);
            mov_rel_sep(j,i,:,:,:) = (squeeze((DataMat(iRep, iTrial, :,:,:))) - img_base)...
                ./img_base;
        end
        
    end
    
%     mov_rel = squeeze(mean(mov_rel_sep, 1)); % do not use this, it squeezes every singleton dimensions (including when trial number = 1)
    mov_rel = mean(mov_rel_sep, 1); mov_rel_size = size(mov_rel);
    mov_rel = reshape(mov_rel, [mov_rel_size(2:end), 1]);
    for i = 1:nPanels
        img_rel(i,:,:) = squeeze(mean(mov_rel(i,:,:,floor(para.fr*opt.tWindow(iTrial, 1))+1 : floor(para.fr*opt.tWindow(iTrial, 2))),4, 'omitnan'));    
    end% data matrix X [#Stim, #Pixel]
    
    X = reshape(img_rel, nPanels, para.width*para.height);
    
    fnametemp = para.filename;
    if ischar(fnametemp)
        vnametemp = [para.pathname{1}, para.filename(1:end-4), '_trial', ...
        num2str(opt.trials(1)),'-',num2str(opt.trials(end)),'_sat=', num2str(opt.ampLimit(2)), '.avi']; 
    else
        vnametemp = [para.pathname{1}, para.filename{1}(1:end-4), '_trial', ...
        num2str(opt.trials(1)),'-',num2str(opt.trials(end)),'_sat=', num2str(opt.ampLimit(2)), '.avi']; 
    end
%     vnametemp = [para.pathname, para.sessionname, '_trial', ...
%         num2str(opt.trials(1)),'-',num2str(opt.trials(end)),'.avi'];

    
%% Plot each repetitions (opt.trials can only be one number)
elseif strcmp(opt.mode, 'allrep')
    nPanels = length(opt.reps) + 1;
    if length(opt.trials) == 1
        iTrial = opt.trials; 
    else
        disp('Error: number of trial selected must be 1')
    end
    
    % ==== calculate relative amplitude change ====
    mov = squeeze(DataMat(opt.reps,iTrial,:,:,:)); 
    mov_rel = mov; % initialize, #Rep x height x width x frames
    img_rel = zeros(length(opt.reps)+1, para.height, para.width); % response mean image for each rep; the last one is averaged
    for i = 1:length(opt.reps)
        mov_temp = squeeze(mov(i,:,:,:)); % 75 x 120 x 25
        img_base = squeeze(mean(mov_temp(:,:,1:floor(para.fr*para.preStim)),3, 'omitnan')); % pre-stimulus: baseline
        img_base = repmat(img_base,1,1,para.nFrame);
        mov_rel(i,:,:,:) = (mov_temp - img_base)./img_base;
        img_rel(i,:,:) = mean(mov_rel(i,:,:,floor(para.fr*opt.tWindow(iTrial, 1))+1 : floor(para.fr*opt.tWindow(iTrial, 2))),4, 'omitnan');
    end
    
    % ==== calculate mean movie, put it after all reps ==== 
    mov_mean = squeeze(mean(mov, 'omitnan')); %averaged across reps, height x width x frames
    img_base = squeeze(mean(mov_mean(:,:,1:floor(para.fr*para.preStim)),3, 'omitnan'));
    img_base = repmat(img_base,1,1,para.nFrame);
    mov_rel(i+1,:,:,:) = (mov_mean - img_base)./img_base;
    img_rel(i+1,:,:) = mean(mov_rel(i+1,:,:,floor(para.fr*opt.tWindow(iTrial, 1))+1 : floor(para.fr*opt.tWindow(iTrial, 2))),4, 'omitnan');
    X = reshape(img_rel, nPanels, para.width*para.height);
    mov_rel_sep = mov_rel;
   
    fnametemp = para.filename;
    if ischar(fnametemp)
        vnametemp = [para.pathname, para.filename(1:end-4), '_trial', ...
        num2str(opt.trials(1)),'-',num2str(opt.trials(end)),'_sat=', num2str(opt.ampLimit(2)), '.avi']; 
    else
        vnametemp = [para.pathname, para.filename{1}(1:end-4), '_trial', ...
        num2str(opt.trials(1)),'-',num2str(opt.trials(end)),'_sat=', num2str(opt.ampLimit(2)), '.avi']; 
    end
end


%% display data and save video
switch opt.plotMode
    
%% combine all subplots into one big matrix (for easier video saving)
    case 'combined'
    if opt.saveON % video-save ON, start a video writer object
        if opt.soundON % video with sound
            
            % ==== for automatic selection, change folder here ==== 
%             folder_sound = 'D:\=code=\McdermottLab\sound_natural\modified_5s\';
%             list = dir(fullfile(folder_sound,'*.wav'));
%             names_sound = natsortfiles({list.name})';
%             [sounddata,fs] = audioread(fullfile(folder_sound,names_sound{opt.trials})); 
            if ~isfield(opt, 'sounddata')
            % ==== for manual selection ====
            [file,path] = uigetfile('*.wav','Select the sound file','D:\SynologyDrive\=sounds=\');
            [sounddata,fs] = audioread(fullfile(path,file)); 
            else 
                fs = 1e5;
                tstart = (para.durStim+para.preStim+para.postStim) * (opt.trials - 1) * fs + 1;
                tend = (para.durStim+para.preStim+para.postStim) * opt.trials * fs;
                sounddata = opt.sounddata(floor(tstart):floor(tend));
            end
            % start a video writer object
            videoFWriter = vision.VideoFileWriter(vnametemp,...
            'FileFormat',       'AVI',...
            'AudioInputPort',   true,...
            'FrameRate',        para.fr,...
            'VideoCompressor',	'None (uncompressed)');   
%             'AudioDataType',    'int16'
            SoundBatchSampleNum =   round(fs/para.fr);
            SoundSeq =              sounddata;
        else % video without sound
            % start a video writer object
            videoFWriter = vision.VideoFileWriter(vnametemp,...
            'FileFormat',       'AVI',...
            'AudioInputPort',   false,...
            'FrameRate',        para.fr,...
            'VideoCompressor',	'None (uncompressed)');
        end
    end

    % arrange subplots 
    if ~isfield(opt, 'p') % automatic calculate #rows and #cols
        [p,n] = numSubplots(nPanels); 
    else % self-defined #rows and #cols
        p = opt.p;
        n = p(1)*p(2);
    end
    
    % construct a combined matrix for mean images
    img_rel = permute(img_rel,[2,3,1]); % height x width x trials
    if nPanels < n % pad with 0 if blocks are left empty
        img_rel(:,:,nPanels+1:n) = zeros(para.height,para.width,n-nPanels);
    end
    
    % rearrange the matrix above into a large matrix
    img_all = [];    
    for k = 1:p(1)
        img_all = [img_all;reshape(img_rel(:,:,p(2)*(k-1)+1:p(2)*k),para.height,para.width*p(2))];
    end

    % FIRST FRAME
    % display the averaged reponse first
    h = imagesc(img_all,opt.ampLimit); axis image; axis off; colorbar;
    colormap(cmap);
    
    pause
%     for i = 0.8*para.fr:4*para.fr % for cropped videos
    for i = 1:size(mov_rel,4)
        mov_temp = mov_rel(:,:,:,i);
        % construct a combined matrix for each frame
        mov_temp = permute(mov_temp,[2,3,1]);
        if nPanels < n % pad with 0 if one block is left empty
            mov_temp(:,:,nPanels+1:n) = zeros(para.height,para.width,n-nPanels);
        end
        mov_all = [];
        for k = 1:p(1)
            mov_all = [mov_all; reshape( mov_temp(:,:,p(2)*(k-1)+1:p(2)*k), para.height, para.width*p(2))];
        end

        set(h,'CData',mov_all)
        title(['time = ',num2str(i*(1./para.fr),'%-5.1f')])
        pause(1/para.fr)

         if opt.saveON   
            if strcmp(opt.color, 'green')
            % save video (in green/black color) and audio
            MAX =           opt.ampLimit(2);
            MIN =           0; % does not show values below zero
            mov_all =       (2^8-1).* (mov_all - MIN)./(MAX - MIN);
            mov_all =       uint8(mov_all);
            mov_all =       repelem(mov_all, 2, 2);% ......... for better video effect
            frame =         repmat(mov_all,1,1,3); 
            frame(:,:,1) =  mov_all*0; % use green color
            frame(:,:,3) =  mov_all*0;
            else 
            % save video (in jet color) and audio
            mov_all = repelem(mov_all, 2, 2);% ......... for better video effect
            Max = opt.ampLimit(2); 
            rgb_ind = floor(mov_all.*(255/2/Max)+257/2);
            rgb_ind(rgb_ind > 256) = 256;
            rgb_ind(rgb_ind < 1) = 1;
            frame_rgb = cmap(rgb_ind,:);
            frame = reshape(frame_rgb, size(mov_all,1), size(mov_all,2), 3);
            end
            
            if opt.soundON
                if i*SoundBatchSampleNum > length(SoundSeq) % if end of current frame is longer than sound (by <1 segment)
                    videoFWriter(frame,SoundSeq((i-1)*SoundBatchSampleNum+1:end) );
                elseif (i-1)*SoundBatchSampleNum+1 > length(SoundSeq) % if sound is already run out
                    videoFWriter(frame);
                else            
                    videoFWriter(frame,SoundSeq((i-1)*SoundBatchSampleNum+1:i*SoundBatchSampleNum) );
                end
            else
                videoFWriter(frame);
            end
        end
    end % end of all frames are processed
    
    % set the last frame as mean response
    set(h,'CData',img_all)
    if opt.saveON % write the last frame and release object videoFWriter
        if strcmp(opt.color, 'green')
        % save video (in green/black color) and audio
        MAX =           opt.ampLimit(2);
        MIN =           0; % does not show values below zero
        img_all =       (2^8-1).*(img_all - MIN)./(MAX - MIN);
        img_all =       uint8(img_all);
        img_all =       repelem(img_all, 2, 2);% ......... for better video effect
%         frame =         ind2rgb(img_all,jet(256));
        frame =         repmat(img_all,1,1,3);
        frame(:,:,1) =  img_all*0;
        frame(:,:,3) =  img_all*0;
          
        else
        % save video (in jet color) and audio
        img_all = repelem(img_all, 2, 2);% ......... for better video effect
        frame_rgb = uint8(zeros(size(img_all,1), size(img_all,2), 3));
        Max = opt.ampLimit(2); 
        rgb_ind = floor(img_all.*(255/2/Max)+257/2);
        rgb_ind(rgb_ind > 256) = 256;
        rgb_ind(rgb_ind < 1) = 1;
        frame_rgb = cmap(rgb_ind,:);
        frame = reshape(frame_rgb, size(img_all,1), size(img_all,2), 3);
        end
        
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
        
        release(videoFWriter) 
    end
    
 %% plot separately, cannot save videos
    case 'separate'
    % show movie as separate matrices (the last is the averaged across reps)
    switch opt.mode
        case 'allrep'
        if ~isfield(opt, 'p')
            [p,n] = numSubplots(nPanels); 
        else
            p = opt.p;
            n = p(1)*p(2);
        end
                
        for iRep = 1:nPanels
            temp = squeeze(mov_rel(iRep,:,:,:));
            subplot(p(1),p(2),iRep); 
            h(iRep) = imagesc(temp(:,:,1),opt.ampLimit); axis image; colorbar; 
            colormap(cmap);
            
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
        for iTrial = 1:nPanels
            temp = squeeze(mov_rel(iTrial,:,:,:));
            subplot(p(1),p(2),iTrial); 
            h(iTrial) = imagesc(temp(:,:,1),opt.ampLimit);colorbar;
            axis image
            title(['Trial.Num.',num2str( opt.trials(iTrial) )])
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