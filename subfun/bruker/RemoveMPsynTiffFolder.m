function [TiffTable,RemovedList] = RemoveMPsynTiffFolder(DataFolder)


TiffFile = dir([DataFolder, '*TSeries*']);
TiffTable=struct2table(TiffFile);
TiffTable=TiffTable(TiffTable.isdir,:);
TiffTable=TiffTable(:,1);
% TiffTable=convertToStrings(TiffTable);

for i=1:size(TiffTable,1)
    tempName=[DataFolder TiffTable.name{i} '\'];
    % [totalRepetitions(i,1), framesAfterStimuli(i,:),laser(i,:),~,~,~,~] = ExpInfoTiffIndiFolder(tempName);
    RemovedList{i}=RemoveMPsynTiff(tempName);
% Zdepth, ZdepthLaser,cycleID,planeID,files,StimID,StimGPLInd
% totalRepetitions, framesAfterStimuli,StimuliPower,Zdepth, ZdepthLaser,cycleID,planeID,files
end
