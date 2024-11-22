
% folderPath='E:\LuSLMOnlineTest\SL0777-Ai203\11142024\Data\TSeries-11142024-0927-028\'
% [totalRepetitions, framesAfterStimuli,StimuliPower,Zdepth, ZdepthLaser,cycleID,planeID,files,StimID,StimGPLInd] = ExpInfoTiffIndiFolder(folderPath)
% 
%%
% 
% RemoveFiles = RemoveMPsynTiff(folderPath)
% 
clear all
ProcessedFolder='E:\LuSLMOnlineTest\SL0777-Ai203\11142024\Data\';
% TiffTable = RemoveMPsynTiffFolder(DataFolder)

[TiffTable, RemoveList] = RemoveMPsynTiffFolder(ProcessedFolder);
