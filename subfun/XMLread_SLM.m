
%Get information from *MarkPoitns.xml files in folder Path
function caTrials=XMLread_SLM(pathGet,caTrialsInd)



%Get information from *MarkPoitns.xml files in folder Path

% pathGet='E:\LuRecording\11022023SLM\';
% [Points,FileInfo,TriggerInfo]=GetMarkPointsInfo(PathF);

caTrialsTemp=dir([pathGet '*TSeries*']);
% caTrialsInd=[5 6];
clear caTrials
for iic = 1:numel(caTrialsInd)
    caTrials(iic) = caTrialsTemp(endsWith({caTrialsTemp.name}, sprintfc('%0.3d',caTrialsInd(iic))));
end




% 
% for iic=1:length(caTrialsInd)
%     iCa=caTrialsInd(icc);
% end


for iic=1:length(caTrialsInd)
    % iic
    clear SeqFile infoAll;
infoAll = xml2struct([caTrials(iic).folder '\' caTrials(iic).name '\' caTrials(iic).name '.xml']);
SeqFile=infoAll.PVScan.Sequence;


clear VolFileInfo MPFileInfo FrameTS LastFame MPrecord numFrames
MPFileInfo=struct([]);
for i=1:length(SeqFile)
    if iscell(SeqFile)
    Temp1=SeqFile{i};
    elseif isstruct(SeqFile)
    Temp1=SeqFile(i);
    else
    end
    if isfield(Temp1,'VoltageRecording')
       % struct2table(Temp1.VoltageRecording.Attributes);
       VolFileInfo(i)=Temp1.VoltageRecording.Attributes;
    end
    if isfield(Temp1,'MarkPoints')
       % struct2table(Temp1.VoltageRecording.Attributes);
       % TempTable1=struct2table(Temp1.MarkPoints.Attributes);
       TempTable1=Temp1.MarkPoints.Attributes;
       if ~isempty(TempTable1.filename)
         
       end

       MPrecord(i)=1;
       if  exist('MPFileInfo')
           % MPFileInfo= concatenateStructs(MPFileInfo,TempTable1);
           % MPFileInfo= [MPFileInfo,TempTable1];
            % MPFileInfo = catstruct(MPFileInfo,TempTable1)
            MPFileInfo = concatenateStructs(MPFileInfo,TempTable1);
       else
           MPFileInfo=TempTable1;
       end
    else
        MPrecord(i)=0;
        % MPFileInfo=[];
    end
    if isfield(Temp1,'Frame')
       clear TempframeTS
       numFrames(i)=length(Temp1.Frame);
       for iS=1:length(Temp1.Frame)
           TempframeTS(iS)=Temp1.Frame{iS}.Attributes;           
       end
       if exist('FrameTS')
          FrameTS=[FrameTS;struct2table(TempframeTS)];
       else
          FrameTS=struct2table(TempframeTS);
       end
    else
       numFrames(i)=0;
    end

end
LastFame=sum(numFrames);
FrameTS = convertTableEntries(FrameTS);
caTrials(iic).FrameTS=FrameTS;
caTrials(iic).FrameNum=numFrames;
caTrials(iic).LastFame=LastFame;
caTrials(iic).MPFileInfo=MPFileInfo;
caTrials(iic).MPrecord=MPrecord;

if exist('VolFileInfo')
caTrials(iic).VolFileInfo=VolFileInfo;



clear vRec
% figure;
for i=1:length(VolFileInfo)
    % if isempty(VolFileInfo(i).dataFile)
    %    break;
    % end

    if isempty(dir(VolFileInfo(i).dataFile))
       break;
   
    end


    vRec{i} = csvread([caTrials(iic).folder '\' caTrials(iic).name '\' VolFileInfo(i).dataFile],1,0);
    numP=size(vRec{i},2)-1;
    % subplot(1,length(SeqFile),i)
    % hold on;
    % for ip=1:numP
    %     plot(vRec{i}(:,1),vRec{i}(:,ip+1));
    % end
end
if exist('vRec')
caTrials(iic).vRec=vRec;
end

clear vConfig
for i=1:length(VolFileInfo)
     if isempty(VolFileInfo(i).configurationFile)
       break
    end
    temp = xml2struct([caTrials(iic).folder '\' caTrials(iic).name '\' VolFileInfo(i).configurationFile]);
    temp1=temp.VRecSessionEntry.Experiment;
    temp2=rmfield(temp.VRecSessionEntry,'Experiment');
    vConfig(i)= MapFields1to2(temp2,temp1);
end
caTrials(iic).vConfig=vConfig;
end


clear pRec TriggerInfo GenInfo pFile AllPoints
PointsInfo=struct([]);
% TriggerInfo=struct([]);
% pRec=struct();
% pFile=struct();
% GenInfo=struct();
for i=1:length(MPFileInfo)
    if isempty(MPFileInfo(i).filename)
       break
    end
    pRec(i) = xml2struct([caTrials(iic).folder '\' caTrials(iic).name '\' MPFileInfo(i).filename]);


    temp=pRec(i).PVMarkPointSeriesElements;
    temp1=temp.Attributes;
    if ~isfield(temp1,'Category')
       temp1.Category=[];
    end
    if ~isfield(temp1,'Name')
       temp1.Name=[];
    end
    GenInfo(i)=temp1;


    temp2=temp.PVMarkPointElement;
    if iscell(temp2)
       temp2=CellTransStruct(temp2);
    end
    for itr=1:length(temp2)
    TriggerInfo{i}{itr}=temp2.Attributes;
    end

    for itr=1:length(temp2)
    tempPoint=temp2(itr).PVGalvoPointElement.Point;
    % pFile(i).Points=struct([]);
    if isstruct(tempPoint)
       tempPoint=StructTransCell(tempPoint);
    end
    for iPoint=1:length(tempPoint)
        if length(PointsInfo)==0
           PointsInfo=tempPoint{iPoint}.Attributes;
        else
           PointsInfo(end+1)=tempPoint{iPoint}.Attributes;
        end
        pFile(i).Points(iPoint)=tempPoint{iPoint}.Attributes;
        % convertTableEntries(struct2table(pFile(1).Points))
    end
    end
    pFile(i).MPtable=convertTableEntries(struct2table(pFile(i).Points));

end
% AllPoints=PointsInfo
% Points=struct2table(PointsInfo);
% uniqueRows = unique(struct2table(PointsInfo), 'rows');
if isempty(MPFileInfo)
   continue
end
AllPoints = convertCellTable(unique(struct2table(PointsInfo), 'rows'));
caTrials(iic).AllPoints=AllPoints;
caTrials(iic).MPTriggerInfo=TriggerInfo;
caTrials(iic).MPRec=pRec;
caTrials(iic).MPFile=pFile;
caTrials(iic).MPGenInfo=GenInfo;




end
