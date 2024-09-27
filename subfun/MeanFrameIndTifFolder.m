function MeanImage=MeanFrameIndTifFolder(TifFolder,nPlane,IndByTotalFrameNum)

% TiffFolder='E:\LuRecording\03082024_MouseMei03_VGlutAi203SLM\TSeries-03082024-1105-003\'
% TiffFile=dir([TiffFolder '*TSeries*Ch2*.ome.tif'])

% % nPlane=3;
% for iPlane=1:nPlane
%     for iframe=1:length(TiffFile)
%         X(:,:,iframe,iPlane)=imread([TiffFile(iframe).folder '\' TiffFile(iframe).name],iPlane);
%     end
%     MeanImage(:,:,iPlane) = mean(X(:,:,:,iPlane),3);
% end
%

% [totalRepetitions, framesAfterStimuli,StimuliPower,Zdepth, ZdepthLaser] = ExpInfoTiffIndiFolder(TifFolder);
[~, ~,~,~, ~,cycleID,planeID,files] = ExpInfoTiffIndiFolder(TifFolder);

if max(planeID)==nPlane
   ComBinedTiff=0;
else 
   ComBinedTiff=1;
   return
end
TotalcycleID=unique(cycleID);
NeedCycleID=TotalcycleID(IndByTotalFrameNum);
% cycleNeed=find(ismember(cycleID,NeedCycleID)==1);
% nPlane=3;
for iframe=1:length(NeedCycleID)
    % TifFile=dir([TifFolder '*TSeries*Cycle* 0' num2str(IndByTotalFrameNum(iframe)) '_Ch2_*.ome.tif']);
    if ComBinedTiff==1
        I1=find(cycleID==NeedCycleID(iframe));
        for iPlane=1:nPlane
            X(:,:,iframe,iPlane)=imread([files(I1).folder '\' files(I1).name],iPlane);
        end
    else
        for iPlane=1:nPlane
            I1=find(cycleID==NeedCycleID(iframe)&planeID==iPlane);
            [files(I1).folder '\' files(I1).name];
            X(:,:,iframe,iPlane)=imread([files(I1).folder '\' files(I1).name],1);
        end
    end
end

for iPlane=1:nPlane
    MeanImage(:,:,iPlane) = mean(X(:,:,:,iPlane),3);
end