function ct = getContour(varargin)
if length(varargin)>=1
    M = varargin{1};
    para = varargin{2};
else
    M = varargin{1};
end

% stimulus parameter for sound "TonePipSeq_A4-A10_73x0.2s(0.2s)@1.0st_(2.7+14.6+2.7)s"
para.tone_dur       = 0.2;% duration of each tone
para.num_per_oct 	= 12; % number of tones per octave
para.tone_gap            = 2.7; % before and after tone pip sequence

mask            = M.saturation>para.sat_lim;
M.hue           = M.hue.*mask;
% map_blur        = imgaussfilt(M.hue,1);

figurex, 
imagesc(M.map_rgb); 
% imagesc(map_hue); 
axis image,

%% get contours, method 1 (using contour function)
% zindex = (para.tone_gap + para.num_per_oct*para.tone_dur: para.num_per_oct*para.tone_dur: para.durStim - para.tone_gap)./...
%     (para.preStim + para.durStim + para.postStim);
zindex = (para.tone_gap : para.num_per_oct*para.tone_dur: para.durStim - para.tone_gap)./...
    (para.preStim + para.durStim + para.postStim);
hold all,
[C,c] = contour(M.hue, zindex, 'LineStyle','-.'); 
axis image; c.LineWidth = 2;     
colormap(hsv), caxis([para.tone_gap para.durStim - para.tone_gap]./para.durStim)

%% generate vectors represent contours
ct = struct;
ct.c = cell(1);
i = 1;
iLevel = 1;
ct.level = [];
while i<= size(C,2)
    ct.level(iLevel) = C(1,i);
    nVer = C(2,i);
    ct.c{iLevel} = C(:,i+1:i+nVer);
    i = i+1+nVer; % the next start point of the sequence
    iLevel = iLevel + 1;   
end
ct.zindex = zindex;
end

