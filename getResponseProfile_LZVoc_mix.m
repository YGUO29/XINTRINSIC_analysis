function [I_inorder, R_inorder, tags_inorder, snames_inorder] = getResponseProfile_LZVoc_mix(R,plot_on)
% filelist = dir('\\10.16.58.229\Test_Procedures\==Slides & Documents\Music\naturalsounds165\naturalsounds165')
folder_origin = 'D:\SynologyDrive\=sounds=\Vocalization\temp_for capsule\AllMix_Orig_norm';
list = dir(fullfile(folder_origin,'*.wav'));
snames = natsortfiles({list.name})';

%% bar plot figure
[R_inorder, I_inorder] = sort(R,'descend');
tags_inorder = R;
snames_inorder = cell(size(tags_inorder));

end