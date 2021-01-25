% =============== load files to analyze ================
% for analyzing colony recording and control sessions
clear all
para.nRep = 0;  
%% read multiple files
[file_mat,path_mat] = uigetfile('*.mat','Select a mat file to analyze','X:\', 'Multiselect', 'on');

if ~iscell(file_mat) % only load one file
        load(fullfile(path_mat,file_mat));
        file_parts = strsplit(file_mat,{'_','.'});
        ind = find(strcmp(file_parts,'P1'));
        nametemp = strjoin(file_parts(1:ind-1),'_');
        load(fullfile(path_mat,[nametemp,'.mat']));
        
        if para.nRep == 0
            para.nRep =     size(P.ProcDataMat,1);
            para.nStim =    size(P.ProcDataMat,2);
            para.height =   size(P.ProcDataMat,3);
            para.width =    size(P.ProcDataMat,4);
            para.nFrame =   size(P.ProcDataMat,5);

            para.fr =       P.ProcFrameRate; % frame rate
            para.preStim =  S.TrlDurPreStim;
            para.durStim =  S.TrlDurStim;
            para.postStim = S.TrlDurPostStim;
%             para.order =    S.SesTrlOrderVec;
            DataMat =       P.ProcDataMat; % DataMat = [rep, trial, height, width, frams]
            para.filename = file_mat;
            para.pathname = path_mat;
        else
            % ==== if need to match speed ====
%             temp = P.ProcDataMat(:,:,:,:,5*para.fr+1:end-5*para.fr);
%             temp = repelem(temp,1,1,1,1,2);
%             DataMat(para.nRep+1 : para.nRep+size(P.ProcDataMat,1),:,:,:,:) = ...
%                 cat(5,P.ProcDataMat(:,:,:,:,1:5*para.fr), temp, P.ProcDataMat(:,:,:,:,end+1-5*para.fr:end));
            % ================================
            DataMat(para.nRep+1 : para.nRep+size(P.ProcDataMat,1),:,:,:,:) = squeeze(P.ProcDataMat);
            para.nRep = para.nRep + size(P.ProcDataMat, 1);
        end
        clear S P
else % load multiple files
    for i = 1:length(file_mat)
        load(fullfile(path_mat,file_mat{i}));
        file_parts = strsplit(file_mat{i},{'_','.'});
        ind = find(strcmp(file_parts,'P1'));
        nametemp = strjoin(file_parts(1:ind-1),'_');
        load(fullfile(path_mat,[nametemp,'.mat']));

    if para.nRep + i == 1 % read the initial file
    %     para.nRep =     size(P.ProcDataMat,1);
        para.nStim =    size(P.ProcDataMat,2);
        para.height =   size(P.ProcDataMat,3);
        para.width =    size(P.ProcDataMat,4);
        para.nFrame =   size(P.ProcDataMat,5);

        para.fr =       P.ProcFrameRate; % frame rate
        para.preStim =  S.TrlDurPreStim;
        para.durStim =  S.TrlDurStim;
        para.postStim = S.TrlDurPostStim;
        para.order(para.nRep+i,:) =    S.SesTrlOrderVec;
        DataMat =       P.ProcDataMat; % DataMat = [rep, trial, height, width, frams]
        para.filename = file_mat;
        para.pathname = path_mat;
    else % continuous reading 
        DataMat(para.nRep+1 : para.nRep+size(P.ProcDataMat,1),:,:,:,:) = squeeze(P.ProcDataMat);
    end
    
    clear S P
    para.nRep = size(DataMat,1);
    disp([num2str(para.nRep),' finished'])
    end     
end

%% get X, and get subtracted DataMat_norm
para.nRep = 20; para.nStim = 2;
DataMat = reshape(DataMat, para.nRep, para.nStim, para.height, para.width, para.nFrame);
% ==== if need to be reversed ====
temp = flip(DataMat(:,2,:,:,:), 5);
DataMat(:,2,:,:,:) = temp;


opt = struct;
[X, DataMat_norm] = getX(DataMat, para, opt);

opt.ampLimit    = 0.4.*[0 1];
% figure, set(gcf, 'color','w')
% [X, DataMat_norm] = ViewData(DataMat, para, opt); % X may con
% DataMat_norm_sub = squeeze(DataMat_norm(1,:,:,:) - DataMat_norm(2,:,:,:));
DataMat_norm_sub = squeeze(DataMat_norm(1,:,:,3:end) - DataMat_norm(2,:,:,1:end-2));

video_sub = DataMat_norm_sub;
fr = para.fr;
%% view video
session_name = 'Reversed';
Max = max(abs(DataMat_norm_sub(:)));
figure, 
imagesc(DataMat_norm_sub(:,:,1), [-Max Max]), axis image; 
title([session_name, ' - Original'])
colorbar
pause
for i = 1:350
    imagesc('CData',DataMat_norm_sub(:,:,i), [-Max Max]), 
    colormap(jet), 
    pause(0.1); 
end
%% plot averaged subtractive image
session_name = '(NR+HP+PS+Rev)';

DataMat_mean = squeeze(mean(DataMat_norm_sub,3));
Max = max(abs(DataMat_norm(:)));

figure, set(gcf, 'color','white','position',[1375 192 560 1075]);
subplot(3,1,1)
imagesc(squeeze(mean(DataMat_norm(1,:,:,:),4)), [-Max, Max]), axis image, colormap(jet), colorbar
title('(NR+HP+PS)')

subplot(3,1,2),
imagesc(squeeze(mean(DataMat_norm(2,:,:,:),4)), [-Max, Max]), axis image, colormap(jet), colorbar
title(session_name)

Max = max(abs(DataMat_mean(:)));
subplot(3,1,3), 
imagesc(DataMat_mean, [-Max, Max]), axis image, colormap(jet), colorbar
title(['(NR+HP+PS) - ', session_name])
set(findall(gcf,'-property','FontSize'),'FontSize',12);

%% read 2 videos, combine them (separately in red and green channels)
filename = cell(1,2);
video = filename;
filename{1} = 'X:\2018.03 T2 (Marmoset 80Z, Xintrinsic, Green & Fluo)\M80Z-180411-12\180411T150747_Blue_Koehler_Fluo_GFP_ColonyOriginal.mp4';
filename{2} = 'X:\2018.03 T2 (Marmoset 80Z, Xintrinsic, Green & Fluo)\M80Z-180411-12\180411T135813_Blue_Koehler_Fluo_GFP_TempoIncreased.mp4';
% read video frames
for i = 1:2
    v = VideoReader(filename{i});
    k = 1;
    while hasFrame(v)
    video{i}(:,:,:,k) = readFrame(v);
    k = k+1;
    end
end
% tempo increased session
% temp = video{2}(:,:,:,5*fr+1:end-5*fr);
% temp = repelem(temp,1,1,1,2);
% video{2} = cat(4,video{2}(:,:,:,1:5*fr), temp, video{2}(:,:,:,end+1-5*fr:end));

video_sub = squeeze(mean(video{2},3)) - squeeze(mean(video{1},3));
fr = v.FrameRate;

% process video frames (put two videos in two color channels)
% video_comb = zeros(size(video{1},1), size(video{1},2), 1400, 3);
% for i = 1:2
%     video_temp = squeeze(mean(video{i},3));
% %     video_temp = repelem(video_temp, 2, 2);
%     video_comb(:,:,:,i) = video_temp;
% end


%% save video (change vname!!)
vname = 'D:\=data=\Marmoset_imaging\video_wf\Vocalization\(NoiseReduce+HighPass+PitchShift)-(NoiseReduce+HighPass+PitchShift+Reversed)_5fps.avi';

% ======= convert to frames (uint8) ======
c = jet(256);
c = uint8(c.*(2^8-1));
video_rgb = uint8(zeros(size(video_sub,1), size(video_sub,2), 3));
video_comb = uint8(zeros(size(video_sub,1), size(video_sub,2), 3, size(video_sub,3)));
Max = max(abs(video_sub(:)));
for i = 1:size(video_comb,4)
video_temp = video_sub(:,:,i);
rgb_ind = floor(video_temp.*(255/2/Max)+257/2);
video_rgb = c(rgb_ind,:);
video_comb(:,:,:,i) = reshape(video_rgb, size(video_sub,1), size(video_sub,2), 3);
end

filename{1} = 'X:\2018.03 T2 (Marmoset 80Z, Xintrinsic, Green & Fluo)\M80Z-180411-12\180411T150747_Blue_Koehler_Fluo_GFP_ColonyOriginal.mp4';
[sounddata, fs] = audioread(filename{1});
videoFWriter = vision.VideoFileWriter(vname,...
            'FileFormat',       'AVI',...
            'AudioInputPort',   true,...
            'FrameRate',        fr,...
            'VideoCompressor',	'None (uncompressed)',...
            'AudioDataType',    'int16');
            SoundBatchSampleNum =   round(fs/fr);
            SoundSeq =              sounddata(:,1);
            
  % save video and audio
% MAX =           0.95* max(video_comb(:));
% MIN =           min(video_comb(:));
% video_comb =       (2^8-1).* (video_comb - MIN)./(MAX - MIN);
% video_comb =       uint8(video_comb);

fwait = waitbar(0,'Start saving frames...');
for i = 1:size(video_comb,4)
    waitbar(i/size(video_comb,4),fwait,['Saving frame ',num2str(i),'/',num2str(size(video_comb,4))]);

    frame = squeeze(video_comb(:,:,:,i));

    if i*SoundBatchSampleNum > length(SoundSeq) % if end of current frame is longer than sound (by <1 segment)
        videoFWriter(frame,SoundSeq((i-1)*SoundBatchSampleNum+1:end) );
    elseif (i-1)*SoundBatchSampleNum+1 > length(SoundSeq) % if sound is already run out
        videoFWriter(frame);
    else            
        videoFWriter(frame,SoundSeq((i-1)*SoundBatchSampleNum+1:i*SoundBatchSampleNum) );
    end

end

release(videoFWriter) 
