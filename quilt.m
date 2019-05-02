% this script demonstrates the Sound Quilting Toolbox
%
% uses default values from Nature Neuroscience (2015) quilting paper
%
% 7-16-2016 -- Josh McDermott (jhm@mit.edu)
% Feb-2019 -- modified YG


%% read sound
clear all,clc
if strcmp(getenv('computername'),'LAPTOP-ERQCNU7B')
    addpath('D:\Dropbox\_Yueqi Guo\_research\Imaging\=code=\McdermottLab\toolbox_sound_quilt_v1.3')
else
    addpath('D:\=code=\McdermottLab\toolbox_sound_quilt_v1.3')
end

folder_origin = 'D:\Dropbox\_Yueqi Guo\_research\Imaging\=sound=\Natural_selected\toilet';
folder_origin_short = [folder_origin,'\origin'];
folder_quilt1 = [folder_origin,'\quilt1'];
folder_quilt2 = [folder_origin,'\quilt2'];
% folder_origin = 'D:\=code=\McdermottLab\sound_texture';
% folder_quilt1 = 'D:\=code=\McdermottLab\sound_texture\quilt1';
if ~exist(folder_quilt1,'dir') || ~exist(folder_quilt2,'dir') || ~exist(folder_origin_short,'dir')
    mkdir(folder_quilt1)
    mkdir(folder_quilt2)
    mkdir(folder_origin_short)
end
% folder_quilt2 = 'D:\=code=\McdermottLab\sound_texture\quilt2';


list = dir(fullfile(folder_origin,'*.wav'));
% ind_empty = find([list.bytes] == 0); % delete empty files
% list(ind_empty) = [];
% which_source = input('Which source file to use (1=speech, 2=writing, 3=wind)? ');
% segment_size_in_win = input('Segment size (number of 30 ms windows, e.g. 2)? ');

%%
segment_size_in_win = 1; %self defined, (number of 30 ms windows, e.g. 2)
%values from NN paper:
win_ms = 30;
stim_length_in_sec = 2; % used to be 6 in original paper

P.filt_density = 1; %1 for regular filterbank, 2 for 4x overcomplete
P.N_audio_channel = 30;
P.audio_low_lim_Hz = 20;
P.audio_high_lim_Hz = 22000; % P.audio_high_lim_Hz - the upper cutoff of the highest filter
P.audio_sr = 44100; % used to be 20000
% for i = 5
%     [source_s,sr] = audioread([folder_origin,'\',list(i).name]);
%     amp_thresh = -35;
%     gap_len_thresh = 8;
%     fname = [];
%     output_dir = '.';
%     display_figures = 1;
%     [new_s,gap_contents] = excise_breaths_and_pauses(source_s,sr,win_ms,amp_thresh,gap_len_thresh,fname,output_dir, display_figures);
% end

%% generate short version of original sound
for i = 1:length(list)
    [source_s,sr] = audioread([folder_origin,'\',list(i).name]);
    if sr~=P.audio_sr
    source_s = resample(source_s,P.audio_sr,sr);
    end
    source_short = source_s(1:P.audio_sr*stim_length_in_sec);
    source_short = source_short./max(abs(source_short));
    filename = [folder_origin_short,'\',list(i).name];
    audiowrite(filename,source_short,P.audio_sr);
end
%% generate quilted sound
for segment_size_in_win = [1 2 4 8 16]
for quilt_number = 1:2
close all, clc    
for i = 1:length(list)
% for i = 3
%     source_s = new_s;
    [source_s,sr] = audioread([folder_origin,'\',list(i).name]);

    if sr~=P.audio_sr
    source_s = resample(source_s,P.audio_sr,sr);
    end
    
%     [source_s, ~] = excise_breaths_and_pauses(source_s,P.audio_sr,win_ms,-35,8,[],[],0);

    [quilted_s, final_seg_order, source_seg_changes, quilt_seg_changes, kl] = ...
    generate_quilt(source_s, win_ms, segment_size_in_win, stim_length_in_sec, P);

    % ==== added by Yueqi 4/2/2019 ====
    quilted_s = quilted_s(floor(1:stim_length_in_sec*P.audio_sr));
    quilted_s = quilted_s./max(abs(quilted_s));
    % ================================
    [n1,bins] = hist(source_seg_changes,40);
    [n2,bins] = hist(quilt_seg_changes,bins);
%     figure; hold on
%     plot(bins,n1/sum(n1),'b','LineWidth',2);
%     plot(bins,n2/sum(n2),'r','LineWidth',2);
%     title('Between-Segment Changes in Original and Quilted Signals');
%     ylabel('Proportion');
%     xlabel('Size of change (arbitrary units)');
%     legend('Original Signal','Quilted Signal');
    
    [~,name_temp,~] = fileparts(list(i).name);
    if quilt_number == 1
        name_temp = [name_temp,'_quilt1'];
        filename = [folder_quilt1,'\',name_temp,'_',num2str(floor(segment_size_in_win*win_ms)),'seg.wav'];
    else
        name_temp = [name_temp,'_quilt2'];
        filename = [folder_quilt2,'\',name_temp,'_',num2str(floor(segment_size_in_win*win_ms)),'seg.wav'];
    end
    audiowrite(filename,quilted_s,P.audio_sr);
end

end
end
% if which_source==1
%     [source_s,sr] = audioread('example_speech_source_file.wav');
% elseif which_source==2
%     [source_s,sr] = audioread('example_pencil_writing_source_file.wav');
% elseif which_source==3
%     [source_s,sr] = audioread('example_wind_source_file.wav');
% end
%% play quilted sound
obj = audioplayer(quilted_s,P.audio_sr);
play(obj)
%% play original sound 
obj = audioplayer(source_s,sr);
play(obj)




