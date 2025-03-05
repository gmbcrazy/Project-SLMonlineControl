
% folderPath='E:\LuSLMOnlineTest\SL0777-Ai203\11142024\Data\TSeries-11142024-0927-028\'
% [totalRepetitions, framesAfterStimuli,StimuliPower,Zdepth, ZdepthLaser,cycleID,planeID,files,StimID,StimGPLInd] = ExpInfoTiffIndiFolder(folderPath)
% 
%%
% 
% RemoveFiles = RemoveMPsynTiff(folderPath)
% 
clear all
ProcessedFolder='E:\LuSLMOnlineTest\SL0838-Ai203\01142025\SingleP\Top14SpeedStimEdgeExc\Data\DataGroup\';
% TiffTable = RemoveMPsynTiffFolder(DataFolder)

RemoveFrame=2; %%2 repetitions right together with MarkPoint would be removed.
[TiffTable, RemoveList] = RemoveMPsynTiffFolder(ProcessedFolder,RemoveFrame);
