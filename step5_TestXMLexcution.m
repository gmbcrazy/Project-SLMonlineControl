clear all
% TestFile='TSeries-04222024-0926-040'
WorkingFolder='E:\LuSLMOnlineTest\SL0838-Ai203\01142025\'
% load('C:\Users\User\Project-SLMonlineControl\subfun\Color\colorMapPN3.mat');
% load('C:\Users\zhangl33\Projects\Project-SLMonlineControl\subfun\Color\colorMapPN3.mat');
confSet = ReadYaml([WorkingFolder 'CurrentSLMsetting.yml']);

ProcessFolder=[WorkingFolder 'SingleP\' 'Top14SpeedStimEdgeExc\'];
% % ConfigFolder='C:\Users\User\Project-SLMonlineControl\config\';
% % 
% % SLMsettingFile='SLMsetting.yml';
% % confSet = ReadYaml([ConfigFolder '\' SLMsettingFile]);

step4_MultiZ_SubStep0_LoadData

umPerPixel=mean([yaml.umPerlPixelX yaml.umPerlPixelY]);



PSTHparam.PreSLMCal=15;        %<----------------------------------------------------------------------------------Edit,Frame # before SLM to calculate baseline map
PSTHparam.PostSLMCal=3;        %<----------------------------------------------------------------------------------Edit,Frame # after SLM to calculate responsive map
PSTHparam.YLim=[-50 600];       % %<-----------------------------------------For Suite2p based ROI signal only method
PSTHparam.pTh=0.05;             % %<-----------------------------------------For Suite2p based ROI signal only method
PSTHparam.TestMethod='ranksum'; % %<-----------------------------------------For Suite2p based ROI signal only method
PSTHparam.MPFrameJump=2;


% idRanges=[3;3];
% idRanges=[4;37];
% idRanges=[FileGenerateInfoTemp.FileID;FileGenerateInfoTemp.FileID];   %Automatic update the new File ID to calculate ROIs
idRanges=[3;5];   %Automatic update the new File ID to calculate ROIs

XMLpattern = 'Laser([\d.]+)GPoint\s?(\d+)';
[ExcuteTBL,MatTBL] = XMLmatch_BinMat(DataFolder, confSet, PSTHparam, Pos3Dneed, idRanges, XMLpattern);
NonMatchFileID=[]
FileIDAll=unique(ExcuteTBL.FileID);

for iFile=1:length(FileIDAll)
     IndExc=find(ExcuteTBL.FileID==FileIDAll(iFile));
     IndMat=find(MatTBL.FileID==FileIDAll(iFile));
     
     ExcuteTBL.Point(IndExc);
     MatTBL.Point(IndMat);

     if length(IndExc)==length(IndMat)
        DiffP=sum(abs(ExcuteTBL.Point(IndExc)-MatTBL.Point(IndMat)));
        if DiffP==0
           disp(['FileID ' num2str(FileIDAll(iFile)) ' xml file excuted correctly']);
        else
           disp(['FileID ' num2str(FileIDAll(iFile)) ' xml file excuted not correctly']);
           disp([MatTBL.Point(IndMat)';ExcuteTBL.Point(IndExc)']);
        end
     else
           disp(['FileID ' num2str(FileIDAll(iFile)) ' bin file < tiff file']);
           disp(MatTBL.Point(IndMat)');
           disp(ExcuteTBL.Point(IndExc)');

     end
end
%% Save results and generate final MarkPoints and Functional groups




SLMPosInfo=load('E:\LuSLMOnlineTest\SL0838-Ai203\01142025\SingleP\Top14SpeedStimEdgeExc\SLMFunGroup.mat');

% SaveFolder='E:\LuSLMOnlineTest\SL0777-Ai203\12182024\SingleP\Top12SpeedStimEdgeExc\';

Pos3DAll=SLMPosInfo.FinalPos3D;
FunScore=SLMPosInfo.FinalFunScore;
Group=SLMPosInfo.Group;
confSet=SLMPosInfo.confSetFinal;
Zdepth=confSet.scan_Z+confSet.ETL

for iGroup=1:length(Group)
    Pos3DGroup{iGroup}=Pos3DAll(Group(iGroup).Indices,:);
end



DataFolder='E:\LuSLMOnlineTest\SL0838-Ai203\01142025\SingleP\Top14SpeedStimEdgeExc\Data\';
idRanges=[38;52];   %Automatic update the new File ID to calculate ROIs
XMLpattern = 'Laser([\d.]+)FunGroup\s?(\d+)';
[ExcuteTBL,MatTBL] = XMLmatch_BinMat(DataFolder, confSet, PSTHparam, Pos3DGroup, idRanges, XMLpattern);
NonMatchFileID=[]
FileIDAll=unique(ExcuteTBL.FileID);

for iFile=1:length(FileIDAll)
     IndExc=find(ExcuteTBL.FileID==FileIDAll(iFile));
     IndMat=find(MatTBL.FileID==FileIDAll(iFile));
     
     ExcuteTBL.Group(IndExc);
     MatTBL.Group(IndMat);

     if length(IndExc)==length(IndMat)
        DiffP=sum(abs(ExcuteTBL.Group(IndExc)-MatTBL.Group(IndMat)));
        if DiffP==0
           disp(['FileID ' num2str(FileIDAll(iFile)) ' xml file excuted correctly']);
        else
           disp(['FileID ' num2str(FileIDAll(iFile)) ' xml file excuted not correctly']);
           disp([MatTBL.Group(IndMat)';ExcuteTBL.Group(IndExc)']);
        end
     else
           disp(['FileID ' num2str(FileIDAll(iFile)) ' bin file < tiff file']);
           disp(MatTBL.Group(IndMat)');
           disp(ExcuteTBL.Group(IndExc)');

     end
end


