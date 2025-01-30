clear all
% TestFile='TSeries-04222024-0926-040'
WorkingFolder='E:\LuSLMOnlineTest\SL0838-Ai203\01292025\'
% load('C:\Users\User\Project-SLMonlineControl\subfun\Color\colorMapPN3.mat');
% load('C:\Users\zhangl33\Projects\Project-SLMonlineControl\subfun\Color\colorMapPN3.mat');
confSet = ReadYaml([WorkingFolder 'CurrentSLMsetting.yml']);

ProcessFolder=[WorkingFolder 'SingleP\' 'Top51SpeedStimEdgeExc\'];
load([ProcessFolder,'SLMFunGroup.mat'])
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
idRanges=[28 36;34 40];   %Automatic update the new File ID to calculate ROIs

XMLpattern = 'Laser([\d.]+)GPoint\s?(\d+)';
[ExcuteTBL,~,MatTBL] = XMLmatch_BinMat(DataFolder, confSet, PSTHparam, Pos3Dneed, idRanges, XMLpattern);
NonMatchFileID=[];

MatTBL.PVpower=round(xmlPower2PVpower(MatTBL.UncagingLaserPower));

FileIDAll=unique(ExcuteTBL.FileID);

a=[ExcuteTBL.FileID ExcuteTBL.UncagingLaserPower MatTBL.PVpower] 

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




% SLMPosInfo=load('E:\LuSLMOnlineTest\SL0838-Ai203\01142025\SingleP\Top14SpeedStimEdgeExc\SLMFunGroup.mat');
SLMPosInfo=load([ProcessFolder,'SLMFunGroup.mat'])

% SaveFolder='E:\LuSLMOnlineTest\SL0777-Ai203\12182024\SingleP\Top12SpeedStimEdgeExc\';

Pos3DAll=SLMPosInfo.FinalPos3D;
FunScore=SLMPosInfo.FinalFunScore;
Group=SLMPosInfo.Group;
confSet=SLMPosInfo.confSetFinal;
Zdepth=confSet.scan_Z+confSet.ETL

for iGroup=1:length(Group)
    Pos3DGroup{iGroup}=Pos3DAll(Group(iGroup).Indices,:);
    Pos3DWeightGroup(:,iGroup)=SLMPosInfo.Group(iGroup).PowerWeight;
end



DataFolder=[ProcessFolder 'Data\'];
idRanges=[42 50;48 76];   %Automatic update the new File ID to calculate ROIs
XMLpattern = 'Laser([\d.]+)FunGroup\s?(\d+)';
[ExcuteTBL, ExcutePowerWeight,MatTBL] = XMLmatch_BinMat(DataFolder, confSet, PSTHparam, Pos3DGroup, idRanges, XMLpattern);

ExcuteWeightGroup=[];
for i=1:size(ExcutePowerWeight)
    minI=find(sum(abs(repmat(ExcutePowerWeight(i,:),length(Group),1)-Pos3DWeightGroup'),2)<0.01);
    if ~isempty(minI)
       ExcuteWeightGroup(i,1)=minI;
    else
       ExcuteWeightGroup(i,1)=0;
    end

end

MatTBL.PVpower=round(xmlPower2PVpower(MatTBL.UncagingLaserPower))

a=[ExcuteTBL.Group ExcuteWeightGroup]


NonMatchFileID=[]
FileIDAll=unique(ExcuteTBL.FileID);

for iFile=1:length(FileIDAll)
     IndExc=find(ExcuteTBL.FileID==FileIDAll(iFile));
     IndMat=find(MatTBL.FileID==FileIDAll(iFile));
     
     ExcuteTBL.Group(IndExc);
     MatTBL.Group(IndMat);

     if length(IndExc)==length(IndMat)
        DiffP=sum(abs(ExcuteTBL.Group(IndExc)-MatTBL.Group(IndMat)));
        DiffPower=sum(abs(ExcuteTBL.UncagingLaserPower(IndExc)-MatTBL.PVpower(IndMat)));
    
        if DiffP==0
           disp(['FileID ' num2str(FileIDAll(iFile)) ' xml file excuted correctly']);
        else
           disp(['FileID ' num2str(FileIDAll(iFile)) ' xml file excuted not correctly']);
           disp([MatTBL.Group(IndMat)';ExcuteTBL.Group(IndExc)']);
        end

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


confSetTemp=confSet;
confSetTemp.save_path0=DataFolder;
[MatTBLvol,FrameTS] = XMLmatch_VoltageOutput(DataFolder, confSetTemp, PSTHparam, idRanges);
disp([num2str(sum(abs(MatTBLvol.VolOut-MatTBLvol.VolOutExc)))  ' errors in WhiskerPulse excution']);


VideoPath=['\\nimhlabstore1.nimh.nih.gov\UFNC\FNC3Z\Camera\SL0838-Ai203\01212025\'];
VideoFile=['16-13-27.seq'];
[roi,roi_mask]=drawBodyROIfromSeq([VideoPath VideoFile],[DataFolder 'Video037']);

DownSampleRByFrame=1;
fps=60;
ShowProcess=0;
fps=100;

fps=83.1315;

[infoTemp, framesWhiskTemp] = Norpix2MATLABAna([VideoPath VideoFile], 10, []);

[diff_hist1,time1]=BodyMovement([VideoPath VideoFile],roi_mask,fps,DownSampleRByFrame,ShowProcess);
[diff_hist1,time1]=BodyMovement([VideoPath VideoFile],roi_mask,fps,DownSampleRByFrame,ShowProcess);

timeX=[1:length(diff_hist1)]/fps;

[fSpeedAll,fStimAll,timeStampCa_PlaneAll,FrameTS,fVideoAll,VideoStartFrameTime]=PV_VolExtract_MultiCyc(confSetTemp,38);

SampleRate=1/median(diff(timeStampCa_PlaneAll(:,1)));


FrameTS.absoluteTime(1)
figure;
plot(FrameTS.absoluteTime(1:3:end)-FrameTS.absoluteTime(1),fSpeedAll(:,1));
hold on;
plot(FrameTS.absoluteTime(1:3:end)-FrameTS.absoluteTime(1),fStimAll(:,1));
% plot(FrameTS.absoluteTime(1:3:end),fVideoAll(:,1)*20);
plot(timeX+VideoStartFrameTime(1),zscore(smooth(diff_hist1)))




figure;
plot(time,diff_hist)
xlabel('Time (Seconds)') 
ylabel('Pixel Change (Intensity)') 


[diff_hist,time]=BodyMovementV0(fileName,roi_mask,fps,DownSampleRByFrame,ShowProcess);









