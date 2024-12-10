%%PostTest Step
clear OutTBL
OutTBL=MPSeqFolder_1TargetXNon(FileGenerateInfo(iCount).tifFolder,confSet,Pos3Dneed);
OutTBL.FileID=FileGenerateInfo(iCount).FileID*ones(size(OutTBL,1),1);
OutTBLAll=[OutTBLAll;OutTBL];
planePos=round(Pos3Dneed(OutTBL.Point,3));
planeDepth=round(confSet.ETL+confSet.scan_Z(1));
[~,planeI]=ismember(planePos,planeDepth);
MPFrameJump=1;
[indexVector, stimulusIDVector, prePostStimuliVector] = getPSTHFrames_MPxmlInfo(OutTBL, PSTHparam.PreSLMCal, PSTHparam.PostSLMCal,MPFrameJump);
PSTHmap = CalMultiPSTHBin(FileGenerateInfo(iCount).binFile, confSet, indexVector, stimulusIDVector, prePostStimuliVector);
PSTHall=cat(3,PSTHall,PSTHmap);

for iP=1:size(OutTBL,1)
    Point=OutTBL.Point(iP);
    roiHeat=Center2neighbour(PSTHmap(:,:,iP,planeI(iP)),Pos3Dneed(Point,[1,2]),ROIparam.NeighbourHfWidthPixel);
    roiHeat=SmoothDecDim3(roiHeat,[1 1])';
       % roiHeat=roiHeat';
    ROIall=cat(3,ROIall,roiHeat);
end
clear PSTHmap roiHeat OutTBL;


FileIDrange=[2;200];                %<------------All files would be used.              
SLMTrialInfo=OutTBLAll;
SLMTrialInfo.Laser=PVpower2xmlPower(SLMTrialInfo.UncagingLaserPower);
[SLMRes,sampleN]=SLMResponse_ROIMultiZ(ROIall,SLMTrialInfo,ROIparam,minTrialN,SumDataFolder,FileIDrange);
iCount=iCount+1;
close all

