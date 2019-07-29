% for comparative analysis
clear all

SessionID = 1; % 1: 80Z; 2: 132D session11; 3: 132D session2; 4: 132D 1st Scrambled; 6: 132D 2nd Scrambled

switch SessionID % 80Z
    case 1
    file_mat = '180724T140401_Blue_Koehler_Fluo_GFP_P1.mat';
    file_tif = '180724T133546_Blue_Koehler_Fluo_GFP.tif';
    path_mat = '\\FANTASIA-DS3617\Test_Imaging\2018.03 T2 (Marmoset 80Z, Xintrinsic, Green & Fluo)\M80Z-180724-13\';
    load(fullfile(path_mat,file_mat));
    [~,nametemp,~,] = fileparts(file_mat);
    load(fullfile(path_mat,[nametemp(1:35),'.mat']));
    I = imread(fullfile(path_mat,file_tif));
    
end

%%