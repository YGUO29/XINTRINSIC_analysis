function [I_inorder, R_inorder, tags_inorder, snames_inorder] = getResponseProfile_NatVoc(R,plot_on)
% filelist = dir('\\10.16.58.229\Test_Procedures\==Slides & Documents\Music\naturalsounds165\naturalsounds165')
folder_origin = 'D:\SynologyDrive\=sounds=\Natural sound\Natural_JM_XINTRINSIC_withLZVoc_200909\Norm';
list = dir(fullfile(folder_origin,'*.wav'));
snames = natsortfiles({list.name})';

% Load Sam's catagory labels directly
load('D:\SynologyDrive\=data=\category_regressors_withLZvoc.mat')
C = C_voc;
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
[R_inorder, I_inorder] = sort(R,'descend');
tags_inorder = R;
snames_inorder = cell(size(tags_inorder));
% close all
if plot_on; 
    figurex; 
end
% ind = [3 2 4 6 5 1]; % for 80Z
% ind = [4 2 3 5 1 6]; % for 132D, session 2
ind = 1:size(R,2);
for i = 1:length(ind)
    if size(R,2) >= ind(i)
        resp = R_inorder(:,ind(i)); index = I_inorder(:,ind(i));
        tags_inorder(:,ind(i)) = tags(index);
        snames_inorder(:,ind(i)) = snames(index);
        Color_inorder = Color(tags_inorder(:,ind(i)),:);
        if plot_on
%             subplot(3,1,i)
            subplot(1,length(ind),i)
            b = bar(resp,'FaceColor','flat');
            title(['Component ',num2str(i)],'fontsize',16)
            b.CData = Color_inorder;
            ylim([-0.5 2.5]), xlim([1 size(R,1)])
            xticks([1 size(R,1)])
            set(gca,'fontsize',24)
            axis square

        end
    else break
    end
end

% if plot_on; figure; end
% for i = 4:6
%     if size(R,2) >= i
%         resp = R_inorder(:,i); index = I_inorder(:,i);
%         tags_inorder(:,i) = tags(index);
%         snames_inorder(:,i) = snames(index);
%         Color_inorder = Color(tags_inorder(:,i),:);
%         if plot_on
% %             subplot(3,1,i-3)
%             subplot(1,6,i)
%             b = bar(resp,'FaceColor','flat');
%             title(['Response Magnitude, Component ',num2str(i)],'fontsize',16)
%             b.CData = Color_inorder;
%             ylim([-0.5 2.5])
%         end
%     else break
%     end
% end

%% plot response amplitude according to catagories
RespGroupMean = zeros(6,nTags);
RespGroupStd  = RespGroupMean;
for t = 1:nTags
    for i = 1:length(ind)
    %     resp = R_inorder(:,i); index = I_inorder(:,i);
    index = find(tags == t);
    resp_temp = R(index,i);
    RespGroupMean(i,t) = mean(resp_temp);
    RespGroupStd(i,t) = std(resp_temp);
    end
end
if plot_on
    figure,
    for i = 1:length(ind)
    subplot(1,length(ind),i),errorbar(RespGroupMean(i,:),RespGroupStd(i,:))
    end
end

end