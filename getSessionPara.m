function para = getSessionPara(para)

% not used anywhere

addpath(genpath(cd))
% para.nRep = 0;  

%% add more repetitions
% 1. run this section to load any number of files from the same folder
% 2. run this section again to select a different set of files... continue
% to stack different repetitions as in "DataMat"

[file_mat,path_mat] = uigetfile('*.mat','Select a mat file to analyze','U:\', 'Multiselect', 'on');
if contains(file_mat, 'GFP')
    mode = 'GFP';
elseif contains(file_mat, 'PBS')
    modality = 'PBS';
else
    disp('check if filename contains GFP or PBS')
end

if ~iscell(file_mat) % only load one file
        load(fullfile(path_mat,file_mat));
        file_parts = strsplit(file_mat,{'_','.'});
%         ind = find(strcmp(file_parts,'GFP')); % GFP or PBS ? 
        ind = find(strcmp(file_parts, modality));
        nametemp = strjoin(file_parts(1:ind),'_');
        load(fullfile(path_mat,[nametemp,'.mat']));
        
        if para.nRep == 0
            para.nRep =     size(P.ProcDataMat,1);
            para.nStim =    size(P.ProcDataMat,2);
            para.height =   size(P.ProcDataMat,3);
            para.width =    size(P.ProcDataMat,4);
            para.nFrame =   size(P.ProcDataMat,5);

            para.fr =       P.ProcFrameRate; % frame rate
            para.preStim =  S.TrlDurPreStim;
            para.durStim =  S.TrlDurStim;
            para.postStim = S.TrlDurPostStim;
            para.order =    S.SesTrlOrderVec;
            DataMat =       P.ProcDataMat; % DataMat = [rep, trial, height, width, frams]
            para.filename = file_mat;
            para.pathname = path_mat;
        else
%             DataMat(para.nRep+1 : para.nRep+size(P.ProcDataMat,1),:,:,:,:) = squeeze(P.ProcDataMat(:,:,6:95,:,:));
            DataMat(para.nRep+1 : para.nRep+size(P.ProcDataMat,1),:,:,:,:) = squeeze(P.ProcDataMat);
            para.nRep = para.nRep + size(P.ProcDataMat, 1);
        end
%         clear S P
else % load multiple files
    for i = 1:length(file_mat)
        load(fullfile(path_mat,file_mat{i}));
        file_parts = strsplit(file_mat{i},{'_','.'});
%         ind = find(strcmp(file_parts,'P1'));
        ind = find(strcmp(file_parts, modality));
        nametemp = strjoin(file_parts(1:ind),'_');
        load(fullfile(path_mat,[nametemp,'.mat']));

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
            para.order(para.nRep+i,:) =    S.SesTrlOrderVec;
            DataMat =       P.ProcDataMat; % DataMat = [rep, trial, height, width, frams]
            para.filename = file_mat;
            para.pathname = path_mat;
        else % continuous reading 
            DataMat(para.nRep+1 : para.nRep+size(P.ProcDataMat,1),:,:,:,:) = squeeze(P.ProcDataMat);
        end
        
%     clear S P
    para.nRep = size(DataMat,1);
    
    disp([num2str(para.nRep),' finished'])
    end     
end
para.side = S.MkySide;
size(DataMat)
para

end