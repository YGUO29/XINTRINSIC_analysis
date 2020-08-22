function getColorbar(Max,mode)

fig1=figure; set(gcf, 'color','w')
axis off
% colormap(jet(100));
colormap(flipud(gray(6)));
if length(Max) == 1
    a = [-Max Max];
else 
    a = Max;
end
caxis(a);
    if strcmp(mode,'percent')
    h = colorbar('Ticks',a,...
        'location','Southoutside',...
        'position',[0.2 0.2 0.5 0.1],...
        'FontSize',26,...
        'Ticklabels',{['-',num2str(Max, '%.1f')], num2str(Max, '%.1f')});        
%         'Ticklabels',{[num2str(100*a(1), '%.1f'), '%'], [num2str(100*a(2), '%.1f'), '%']});
    elseif strcmp(mode,'number')
        h = colorbar('Ticks',a,...
        'location','Southoutside',...
        'position',[0.2 0.2 0.5 0.1],...
        'FontSize',26,...
        'Ticklabels',{num2str(a(1), '%.1f'), num2str(a(2), '%.1f')});
    else
        disp('set mode as percent or number')
    end
end