function [DataMat, para] = readSession(opt)
para.nRep = 0;

% ========== 102D ============
if strcmp(opt.animal, '102D') && strcmp(opt.session, 'NatVocMM') && ...
        strcmp(opt.modal, 'Calcium') && strcmp(opt.date, '201226')
    
    para.pathname = [
        repmat({'U:\2019.09.T1 (Marmoset 102D, Xintrinsic)\M102D-201226-11'}, 1, 6),...
        repmat({'U:\2019.09.T1 (Marmoset 102D, Xintrinsic)\M102D-201227-12'}, 1, 6)]';
    para.filename = {
        '201226T122449_Blue_Koehler_Fluo_GFP_90x144@5fps_P1.mat',...
        '201226T125511_Blue_Koehler_Fluo_GFP_90x144@5fps_P1.mat',...
        '201226T132658_Blue_Koehler_Fluo_GFP_90x144@5fps_P1.mat',...
        '201226T135806_Blue_Koehler_Fluo_GFP_90x144@5fps_P1.mat',...
        '201226T142952_Blue_Koehler_Fluo_GFP_90x144@5fps_P1.mat',...
        '201226T150014_Blue_Koehler_Fluo_GFP_90x144@5fps_P1.mat',...
        '201227T131246_Blue_Koehler_Fluo_GFP_90x144@5fps_P1.mat',...
        '201227T134333_Blue_Koehler_Fluo_GFP_90x144@5fps_P1.mat',...
        '201227T141426_Blue_Koehler_Fluo_GFP_90x144@5fps_P1.mat',...
        '201227T144508_Blue_Koehler_Fluo_GFP_90x144@5fps_P1.mat',...
        '201227T151639_Blue_Koehler_Fluo_GFP_90x144@5fps_P1.mat',...
        '201227T154756_Blue_Koehler_Fluo_GFP_90x144@5fps_P1.mat'}';
end

if strcmp(opt.animal, '102D') && strcmp(opt.session, 'VocMM') && ...
        strcmp(opt.modal, 'Calcium') && strcmp(opt.date, '210511')
    
    para.pathname = [
        repmat({'U:\2019.09.T1 (Marmoset 102D, Xintrinsic)\M102D-210511-12'}, 1, 2),...
        repmat({'U:\2019.09.T1 (Marmoset 102D, Xintrinsic)\M102D-210513-10'}, 1, 8)]';
    para.filename = {
        '210511T140616_Blue_Koehler_Fluo_GFP_90x144@5fps_P1.mat',...
        '210511T143125_Blue_Koehler_Fluo_GFP_90x144@5fps_P1.mat',...
        '210513T105357_Blue_Koehler_Fluo_GFP_90x144@5fps_P1.mat',...
        '210513T112817_Blue_Koehler_Fluo_GFP_90x144@5fps_P1.mat',...
        '210513T115801_Blue_Koehler_Fluo_GFP_90x144@5fps_P1.mat',...
        '210513T122510_Blue_Koehler_Fluo_GFP_90x144@5fps_P1.mat',...
        '210513T125204_Blue_Koehler_Fluo_GFP_90x144@5fps_P1.mat',...
        '210513T131805_Blue_Koehler_Fluo_GFP_90x144@5fps_P1.mat',...
        '210513T134413_Blue_Koehler_Fluo_GFP_90x144@5fps_P1.mat',...
        '210513T141013_Blue_Koehler_Fluo_GFP_90x144@5fps_P1.mat'}';
end


if strcmp(opt.animal, '102D') && strcmp(opt.session, 'Voc') && ...
        strcmp(opt.modal, 'Calcium') && strcmp(opt.date, '210121')
    
    para.pathname = [repmat({'U:\2019.09.T1 (Marmoset 102D, Xintrinsic)\M102D-210121-09'}, 1, 7)]';
    para.filename = {
        '210121T113754_Blue_Koehler_Fluo_GFP_90x144@20fps_P1.mat',...
        '210121T122650_Blue_Koehler_Fluo_GFP_90x144@20fps_P1.mat',...
        '210121T125425_Blue_Koehler_Fluo_GFP_90x144@20fps_P1.mat',...
        '210121T134331_Blue_Koehler_Fluo_GFP_90x144@20fps_P1.mat',...
        '210121T141035_Blue_Koehler_Fluo_GFP_90x144@20fps_P1.mat',...
        '210121T143926_Blue_Koehler_Fluo_GFP_90x144@20fps_P1.mat',...
        '210121T150641_Blue_Koehler_Fluo_GFP_90x144@20fps_P1.mat'}';
end

if strcmp(opt.animal, '102D') && strcmp(opt.session, 'Ripple') && ...
        strcmp(opt.modal, 'Calcium') && strcmp(opt.date, '201211')
    
    para.pathname = [
        repmat({'U:\2019.09.T1 (Marmoset 102D, Xintrinsic)\M102D-201211-10'}, 1, 8),...
        repmat({'U:\2019.09.T1 (Marmoset 102D, Xintrinsic)\M102D-201216-09'}, 1, 3)]';
    para.filename = {
        '201211T125641_Blue_Koehler_Fluo_GFP_90x144@20fps_P1.mat',...
        '201211T131130_Blue_Koehler_Fluo_GFP_90x144@20fps_P1.mat',...
        '201211T132443_Blue_Koehler_Fluo_GFP_90x144@20fps_P1.mat',...
        '201211T133902_Blue_Koehler_Fluo_GFP_90x144@20fps_P1.mat',...
        '201211T135222_Blue_Koehler_Fluo_GFP_90x144@20fps_P1.mat',...
        '201211T140705_Blue_Koehler_Fluo_GFP_90x144@20fps_P1.mat',...
        '201211T142137_Blue_Koehler_Fluo_GFP_90x144@20fps_P1.mat',...
        '201211T143545_Blue_Koehler_Fluo_GFP_90x144@20fps_P1.mat',...
        '201216T102744_Blue_Koehler_Fluo_GFP_90x144@20fps_P1.mat',...
        '201216T110645_Blue_Koehler_Fluo_GFP_90x144@20fps_P1.mat',...
        '201216T112313_Blue_Koehler_Fluo_GFP_90x144@20fps_P1.mat'}';
end

if strcmp(opt.animal, '102D') && strcmp(opt.session, 'Tone') && ...
        strcmp(opt.modal, 'Calcium') && strcmp(opt.date, '201227')
    para.pathname = {'U:\2019.09.T1 (Marmoset 102D, Xintrinsic)\M102D-201227-12'};
    para.filename = {'201227T161929_Blue_Koehler_Fluo_GFP_90x144@20fps_P1.mat'};
end   

% ========== 80Z ============
if strcmp(opt.animal, '80Z') && strcmp(opt.session, 'Nat') && ...
        strcmp(opt.modal, 'Calcium') && strcmp(opt.date, '180724')
    para.pathname = {'Y:\2018.03_M80Z\M80Z-180724-13'};
    para.filename = {'180724T140401_Blue_Koehler_Fluo_GFP_150x240@5fps_P1.mat'};
end

if strcmp(opt.animal, '80Z') && strcmp(opt.session, 'Tone') && ...
        strcmp(opt.modal, 'Calcium') && strcmp(opt.date, '180726')
    para.pathname = {'Y:\2018.03_M80Z\M80Z-180726-14'};
    para.filename = {'180726T144709_Blue_Koehler_Fluo_GFP_150x240@5fps_P1.mat'};
end

% ========== 132D ============
if strcmp(opt.animal, '132D') && strcmp(opt.session, 'Nat') && ...
        strcmp(opt.modal, 'Calcium') && strcmp(opt.date, '190119')
    
    para.pathname = {'Y:\2018.11_M132D\M132D-190119-10'};
    para.filename = {'190119T114434_Blue_Koehler_Fluo_GFP_150x240@5fps_P1.mat'};
end

if strcmp(opt.animal, '132D') && strcmp(opt.session, 'Nat') && ...
        strcmp(opt.modal, 'Calcium') && strcmp(opt.date, '190224')
    
    para.pathname = {'Y:\2018.11_M132D\M132D-190224-09'};
    para.filename = {'190224T095545_Blue_Koehler_Fluo_GFP_150x240@5fps_P1.mat'};
end

if strcmp(opt.animal, '132D') && strcmp(opt.session, 'Tone') && ...
        strcmp(opt.modal, 'Calcium') && strcmp(opt.date, '190409')
    para.pathname = {'Y:\2018.11_M132D\M132D-190409-10'};
    para.filename = {'190409T115402_Blue_Koehler_Fluo_GFP_150x240@5fps_P1.mat'};
end


% ===== load files =====
for i = 1:length(para.filename)
    load(fullfile(para.pathname{i},para.filename{i}));

    file_parts = strsplit(para.filename{i},{'_','.'});
%         ind = find(strcmp(file_parts,'P1'));
    ind = find(strcmp(file_parts, 'P1'));
    nametemp = strjoin(file_parts(1:ind-2),'_');
    load(fullfile(para.pathname{i},[nametemp,'.mat']));

    if para.nRep + i == 1 % read the initial file, initialize parameters
    %     para.nRep =     size(P.ProcDataMat,1);
        para.nStim =    size(P.ProcDataMat,2);
        para.height =   size(P.ProcDataMat,3);
        para.width =    size(P.ProcDataMat,4);
        para.nFrame =   size(P.ProcDataMat,5);
        para.fr =       P.ProcFrameRate; % frame rate

        para.preStim =  S.TrlDurPreStim;
        para.durStim =  S.TrlDurStim;
        para.postStim = S.TrlDurPostStim;
        DataMat =       P.ProcDataMat; % DataMat = [rep, trial, height, width, frams]
    else % continuous reading 
        DataMat(para.nRep+1 : para.nRep + size(P.ProcDataMat,1),:,:,:,:) = squeeze(P.ProcDataMat);
    end
    
    para.order(para.nRep+1 : para.nRep + size(P.ProcDataMat,1),:) ...
        = S.SesTrlOrderMat;


    para.nRep = size(DataMat,1);
    disp(['Repetition ', num2str(para.nRep),' finished'])
end     

% ==== read subset of the data ====
% for 102D: NatVocMM session
ind_voc_orig = 1:2:30; ind_voc_synth = 2:2:30; 
ind_nat_orig = 31:2:size(DataMat,2); ind_nat_synth = 32:2:size(DataMat,2); 

if isfield(opt,'subset') && strcmp(opt.subset,'Nat') %only read natural sound trials
    DataMat = DataMat(:, ind_nat_orig, :, :, :);
elseif isfield(opt,'subset') && strcmp(opt.subset,'NatMM')% only read natural sound + model matched natural sound
    DataMat = DataMat(:, sort([ind_nat_orig, ind_nat_synth]), :, :, :);
elseif isfield(opt,'subset') && strcmp(opt.subset,'VocMM')% only read vocalizations + model matched vocalizations
    DataMat = DataMat(:, sort([ind_voc_synth, ind_voc_orig]), :, :, :);
elseif isfield(opt,'subset') && strcmp(opt.subset,'NatVoc')% only read natural sounds + vocalizations
    DataMat = DataMat(:, sort([ind_nat_orig, ind_voc_orig]), :, :, :);
else
end

% for 80Z, tone session:
if isfield(opt,'subset') && strcmp(opt.subset,'50dB') %only read natural sound trials
    DataMat = DataMat(:, 3*9+1:3*9+9, :, :, :);
end

para.nStim = size(DataMat,2);
end