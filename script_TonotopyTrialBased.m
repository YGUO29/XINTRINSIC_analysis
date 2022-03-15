%% for function version, see getTonotopy.m
%%
% plot trial-based tonotopy
% display mode: 1 = winner take all; 2 = weighted average
display_mode = 2;
smoothON = 1;
smooth_window = 10;
freq_per_oct = 12;

% For sound 'Tones_A2-A10_97x1st_0.6+0.6+1.8s_50dBSPL_ramped_201227.wav'
octaves = 0:1/freq_per_oct:8; % one semitone apart
freqs = 110.*2.^octaves; % A2 ~ A10

if smoothON
    [max_amp, max_ind] = max(smoothdata(X,1,'movmean',smooth_window), [], 1); % find max amplitude (and frequency index) for each pixel
else
    [max_amp, max_ind] = max(X, [], 1); % find max amplitude (and frequency index) for each pixel
end

% under two display mode: use two different scaling of hsv colorbar
C = cell(1,2);
% C{1} = 0:0.01:0.9; % cut the end a little bit because it overlaps with the red in the beginning
C{2} = HueRedist('HSVadjusted', 'Circular'); % 299 numbers in HSV space
C{1} = C{2}(1:270); % cut the end a little bit because it overlaps with the red in the beginning

switch display_mode 
    case 1
        
        bfs = reshape(max_ind, para.height, para.width);
        cmap = interp1(0:length(C{display_mode})-1, C{display_mode},...
            linspace(0, length(C{display_mode})-1, length(freqs))');
        cmap = hsv2rgb([cmap, ones(size(cmap)), ones(size(cmap))]);
%         cmap = hsv(length(freqs));
        tonotopy_rgb = cmap(bfs, :);
    case 2
        X_rectified = X;
        X_rectified(X<0) = 0;
        weights = X_rectified./sum(X_rectified,1);
        
        weighted_avg = freqs*weights; % units = Hz
        weighted_avg = log2(weighted_avg./110); % units = octave above A2

        % convert BF to a color index
        amp_lim = [0 8]; % units: octave above A2
        ind = (length(C{display_mode}) - 1).*(weighted_avg - amp_lim(1))./(amp_lim(2) - amp_lim(1)) + 1;

        cmap = hsv2rgb([C{display_mode}, ones(size(C{display_mode})), ones(size(C{display_mode}))]);
        ind(ismissing(ind)) = 1;
        
        tonotopy_rgb = cmap(floor(ind), :);
    otherwise
end

% plot mode 2: BFs weighted with max response amplitude
mask_threshold = 0.02;
tonotopy_mask = max_amp > mask_threshold;
tonotopy_mask = repmat(tonotopy_mask, 3, 1)';
tonotopy_masked = tonotopy_rgb .* tonotopy_mask;
tonotopy_masked = (tonotopy_masked - min(tonotopy_masked(:)))./(max(tonotopy_masked(:)) - min(tonotopy_masked(:)));

% plot mode 3: BFs with an amplitude mask
sat_threshold = 0.08;
max_amp(max_amp > sat_threshold) = sat_threshold;
tonotopy_amp = repmat(max_amp, 3, 1)';
tonotopy_final = tonotopy_rgb .* tonotopy_amp;
tonotopy_final = (tonotopy_final - min(tonotopy_final(:)))./(max(tonotopy_final(:)) - min(tonotopy_final(:)));

tonotopy = cell(1,3);
tonotopy{1} = reshape(tonotopy_rgb, para.height, para.width, 3); % only BFs
tonotopy{2} = reshape(tonotopy_masked, para.height, para.width, 3); % BFs with an amplitude mask
tonotopy{3} = reshape(tonotopy_final, para.height, para.width, 3); % BFs weighted with max response amplitude

figurex;
Titles = {'BFs only', ['BFs with amplitude < ', num2str(100*mask_threshold), '% mask'], ...
    ['BFs weighted with max response amplitude (sat=', num2str(100*sat_threshold), '%)']};
for i = 1:3 % 3 plot mode
    subplot(1,3,i)
    img = reshape(tonotopy{i}, para.height, para.width, 3);
    imagesc(img), axis image, axis off
    title(Titles{i})
    colormap(cmap)
    if display_mode == 1
        colorbar('Ticks',[1/length(freqs):freq_per_oct/length(freqs):1],...
         'TickLabels',{'A2', 'A3', 'A4', 'A5', 'A6', 'A7', 'A8', 'A9', 'A10'});
    else
        colorbar('Ticks',[1/8:1/8:1],...
             'TickLabels',{'A3', 'A4', 'A5', 'A6', 'A7', 'A8', 'A9', 'A10'});
    end
     
end

%% change colorbar to start from A4
% plot trial-based tonotopy
% display mode: 1 = winner take all; 2 = weighted average
display_mode = 1;
smoothON = 3;
smooth_window = 3;
freq_per_oct = 12;


% octaves = 0:8; % one octave apart
octaves = 0:1/freq_per_oct:8+1/4; % one semitone apart
freqs = 110.*2.^octaves;
if smoothON
    [max_amp, max_ind] = max(smoothdata(X,1,'movmean',smooth_window), [], 1); % find max amplitude (and frequency index) for each pixel
else
    [max_amp, max_ind] = max(X, [], 1); % find max amplitude (and frequency index) for each pixel
end
% under two display mode: use two different scaling of hsv colorbar
C = cell(1,2);
% C{1} = 0:0.01:0.9; % cut the end a little bit because it overlaps with the red in the beginning
C{2} = HueRedist('HSVadjusted', 'Circular'); % 299 numbers in HSV space
C{1} = C{2}(1:270); % cut the end a little bit because it overlaps with the red in the beginning


switch display_mode 
    case 1
        % cut the lowest 2 octaves out of 8 octaves(to match cycle based tonotopy, which starts
        % from A4) 
        max_ind(max_ind < 2*freq_per_oct + 1) = 2*freq_per_oct + 1;
        max_ind = max_ind - 2*freq_per_oct;
        freqs = freqs(2*freq_per_oct + 1:end);
        
        bfs = reshape(max_ind, para.height, para.width);
%         cmap = hsv(length(freqs));
        cmap = interp1(0:length(C{display_mode})-1, C{display_mode}, linspace(0, length(C{display_mode})-1, length(freqs))');
        cmap = hsv2rgb([cmap, ones(size(cmap)), ones(size(cmap))]);
        tonotopy_rgb = cmap(bfs, :);
    case 2
        X_rectified = X;
        X_rectified(X<0) = 0;
        weights = X_rectified./sum(X_rectified,1);
        
        
        weighted_avg = freqs*weights; % units = Hz
        weighted_avg = log2(weighted_avg./110); % units = octave above A2
        % change colorbar range
        weighted_avg(weighted_avg < 2) = 2;
%         cmap = hsv(256);
        % convert BF to a color index
%         ind = (length(cmap) - 1).*rescale(weighted_avg, 0, 1) + 1;
        amp_lim = [2 8]; % units: octave above A2
        ind = (length(C{display_mode}) - 1).*(weighted_avg - amp_lim(1))./(amp_lim(2) - amp_lim(1)) + 1;
        cmap = hsv2rgb([C{display_mode}, ones(size(C{display_mode})), ones(size(C{display_mode}))]);
        ind(ismissing(ind)) = 1;
        
        tonotopy_rgb = cmap(floor(ind), :);
    otherwise
end

% plot mode 2: BFs weighted with max response amplitude
mask_threshold = 0.02;
tonotopy_mask = max_amp > mask_threshold;
tonotopy_mask = repmat(tonotopy_mask, 3, 1)';
tonotopy_masked = tonotopy_rgb .* tonotopy_mask;
tonotopy_masked = (tonotopy_masked - min(tonotopy_masked(:)))./(max(tonotopy_masked(:)) - min(tonotopy_masked(:)));

% plot mode 3: BFs with an amplitude mask
sat_threshold = 0.05;
max_amp(max_amp > sat_threshold) = sat_threshold;
tonotopy_amp = repmat(max_amp, 3, 1)';
tonotopy_final = tonotopy_rgb .* tonotopy_amp;
tonotopy_final = (tonotopy_final - min(tonotopy_final(:)))./(max(tonotopy_final(:)) - min(tonotopy_final(:)));

tonotopy = cell(1,3);
tonotopy{1} = reshape(tonotopy_rgb, para.height, para.width, 3); % only BFs
tonotopy{2} = reshape(tonotopy_masked, para.height, para.width, 3); % BFs with an amplitude mask
tonotopy{3} = reshape(tonotopy_final, para.height, para.width, 3); % BFs weighted with max response amplitude

figurex;
Titles = {'BFs only', ['BFs with amplitude < ', num2str(100*mask_threshold), '% mask'], ...
    ['BFs weighted with max response amplitude (sat=', num2str(100*sat_threshold), '%)']};
for i = 1:3 % 3 plot mode
    subplot(1,3,i)
    img = reshape(tonotopy{i}, para.height, para.width, 3);
    imagesc(img), axis image
    title(Titles{i})
    colormap(cmap)
%     colorbar('Ticks',[1/8:1/8:1],...
%          'TickLabels',{'A3', 'A4', 'A5', 'A6', 'A7', 'A8', 'A9', 'A10'});
    if display_mode == 1
        colorbar('Ticks',[1/length(freqs):freq_per_oct/length(freqs):1],...
         'TickLabels',{'A4', 'A5', 'A6', 'A7', 'A8', 'A9', 'A10'});
%         colorbar('Ticks',[1/7:1/7:1],...
%          'TickLabels',{'A4', 'A5', 'A6', 'A7', 'A8', 'A9', 'A10'});
    else
        colorbar('Ticks',[1/6:1/6:1],...
         'TickLabels',{'A5', 'A6', 'A7', 'A8', 'A9', 'A10'});
    end
end

%%
save('E:\Dropbox\_Yueqi Guo\_research\Imaging\=data=\80Z\180726_tonotopy_Weighted_A4-A10', 'tonotopy')