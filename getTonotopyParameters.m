function para = getTonotopyParameters(para, animal)
switch animal
    case '102D'
        para.display_mode = 1; %1 = winner take all; 2 = weighted average
        para.smoothON = 1;
        para.smooth_window = 10;
        para.freq_per_oct = 12;
        para.freq_base = 110;
        para.freq_oct_range = 8; % A2~C11: 8+1/4; A2~A10: 8;
        para.freq_cut_oct = 0; % cut lower octaves for display purpose
    case {'80Z', '132D'}
        para.display_mode = 2; %1 = winner take all; 2 = weighted average
        para.smoothON = 0;
        para.smooth_window = 10;
        para.freq_per_oct = 1;
        para.freq_base = 110;
        para.freq_oct_range = 8; % A2~C11: 8+1/4; A2~A10: 8;
        para.freq_cut_oct = 2; % cut lower octaves for display purpose
    otherwise
end
  
