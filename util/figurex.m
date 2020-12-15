function h = figurex(p)
h = figure;
set(gcf,...
    'DefaultAxesFontSize',18,...
    'DefaultLineLineWidth', 2,...
    'color','w');
if nargin ~= 0
    set(gcf,'position',p);
else
end
end