function ct = getContour_compmponent(W, comp)
% W: matrix - #comp x #pix
% comp: 3d matrix - height x width x #comp
% ct: cell of size 1 x #comp
% use plotContour.m to plot it

nComp = size(W,1);
ct = struct;
for iComp = 1:nComp
    temp = W(iComp,:);
    max_amp = max(abs(temp));
    zindex = [-max_amp.*[1:-0.1:0.3], max_amp.*[0.3:0.1:1]];
%     zindex = prctile(abs(temp),20:10:90);
    [C, ~] = contour(comp(:,:,iComp), zindex, 'LineStyle','-.'); 
    % generate vectors represent contours
    ct(iComp).c = cell(1);
    i = 1;
    iLevel = 1;
    ct(iComp).level = [];
    while i<= size(C,2)
        ct(iComp).level(iLevel) = C(1,i); % first colume, first row: level
        nVer = C(2,i); % first colume, second row: number of vertices
        ct(iComp).c{iLevel} = C(:,i+1:i+nVer); % read vertices for this level
        i = i+1+nVer; % the next start point of the sequence
        iLevel = iLevel + 1;   
    end
    ct(iComp).zindex = zindex;
end
end

