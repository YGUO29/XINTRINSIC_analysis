function [I_inorder, R_inorder, tags_inorder, snames_inorder] = getResponseProfile(R,plot_on)
% filelist = dir('\\10.16.58.229\Test_Procedures\==Slides & Documents\Music\naturalsounds165\naturalsounds165')

folder_origin = 'D:\=code=\McdermottLab\sound_natural';
list = dir(fullfile(folder_origin,'*.wav'));
snames = natsortfiles({list.name})';

% Load Sam's catagory labels directly
load('D:\=code=\XINTRINSIC_analysis\category_regressors.mat')
tags = C.category_assignments; 
nTags = max(tags);
Color = C.colors;


% [tags,snames] = xlsread([folder_origin,'\NatSound_label'],1);
% set colors (color samples from Sam's 2015 papper by Yueqi)
% Color(1,:) = [19  78  150]./255;
% Color(2,:) = [0   153 211]./255;
% Color(3,:) = [0   104 78 ]./255;
% Color(4,:) = [48  191 159]./255;
% Color(5,:) = [124 66  150]./255;
% Color(6,:) = [162 121 186]./255;
% Color(7,:) = [202 51  32 ]./255;
% Color(8,:) = [255 103 132]./255;
% Color(9,:) = [119 119 119]./255;
% Color(10,:) = [255 139 24]./255;
% Color(11,:) = [255 255 0]./255;

%% bar plot figure
[R_inorder,I_inorder] = sort(R,'descend');
tags_inorder = R;
snames_inorder = cell(size(tags_inorder));
close all
if plot_on; figure; end
for i = 1:3
    if size(R,2) >= i
        resp = R_inorder(:,i); index = I_inorder(:,i);
        tags_inorder(:,i) = tags(index);
        snames_inorder(:,i) = snames(index);
        Color_inorder = Color(tags_inorder(:,i),:);
        if plot_on
            subplot(3,1,i)
            b = bar(resp,'FaceColor','flat');
            title(['Response Magnitude, Component ',num2str(i)],'fontsize',16)
            b.CData = Color_inorder;
        end
    else break
    end
end

if plot_on; figure; end
for i = 4:6
    if size(R,2) >= i
        resp = R_inorder(:,i); index = I_inorder(:,i);
        tags_inorder(:,i) = tags(index);
        snames_inorder(:,i) = snames(index);
        Color_inorder = Color(tags_inorder(:,i),:);
        if plot_on
            subplot(3,1,i-3)
            b = bar(resp,'FaceColor','flat');
            title(['Response Magnitude, Component ',num2str(i)],'fontsize',16)
            b.CData = Color_inorder;
        end
    else break
    end
end

%% plot response amplitude according to catagories
RespGroupMean = zeros(6,nTags);
RespGroupStd  = RespGroupMean;
for t = 1:nTags
    for i = 1:6
    %     resp = R_inorder(:,i); index = I_inorder(:,i);
    index = find(tags == t);
    resp_temp = R(index,i);
    RespGroupMean(i,t) = mean(resp_temp);
    RespGroupStd(i,t) = std(resp_temp);
    end
end
if plot_on
    figure,
    for i = 1:6
    subplot(2,3,i),errorbar(RespGroupMean(i,:),RespGroupStd(i,:))
    end
end

end