function [DataMat, para] = getDataMat
para.nRep = 0;
[file_mat,path_mat] = uigetfile('*.mat','Select a mat file to analyze','X:\', 'Multiselect', 'on');

if ~iscell(file_mat) % only load one file
        load(fullfile(path_mat,file_mat));
        file_parts = strsplit(file_mat,{'_','.'});
        ind = find(strcmp(file_parts,'P1'));
        nametemp = strjoin(file_parts(1:ind-1),'_');
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
%             para.order =    S.SesTrlOrderVec;
            DataMat =       P.ProcDataMat; % DataMat = [rep, trial, height, width, frams]
            para.filename = file_mat;
            para.pathname = path_mat;
        else
            DataMat(para.nRep+1 : para.nRep+size(P.ProcDataMat,1),:,:,:,:) = squeeze(P.ProcDataMat);
            para.nRep = para.nRep + 1;
        end
        clear S P
else % load multiple files
    for i = 1:length(file_mat)
        load(fullfile(path_mat,file_mat{i}));
        file_parts = strsplit(file_mat{i},{'_','.'});
        ind = find(strcmp(file_parts,'P1'));
        nametemp = strjoin(file_parts(1:ind-1),'_');
        load(fullfile(path_mat,[nametemp,'.mat']));

    if para.nRep + i == 1 % read the initial file
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
    
    clear S P
    para.nRep = size(DataMat,1);
    disp([num2str(para.nRep),' finished'])
    end     
end

end
