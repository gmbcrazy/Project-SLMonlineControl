clear all

WorkFolder='E:\LuSLMOnlineTest\SL0838-Ai203\01292025\';
ConfigFile='CurrentSLMsetting.yml';%<----------------------------------------------------------------------------------Edit, configuration file
[~,~,~,CaData,CaDataPlane,stat,yaml,confSet]=ROIToXYZ(WorkFolder,ConfigFile);
ProcessFolder = Get_ExpDataFolder(WorkFolder, 'SpeedStimEdgeExc', {'Data','AllIncluded','DataSum','.gpl','.xml'})

SLMPosInfo=load([ProcessFolder 'SLMFunGroup.mat']);
load('C:\Users\zhangl33\Projects\Project-SLMonlineControl\CodeAndPackages\subfun\Color\colorMapPN3.mat','colorMapPN1')

DataFolder=[ProcessFolder 'Data\'];

SaveFolder=[ProcessFolder 'Results\'];
mkdir(SaveFolder)

Pos3DAll=SLMPosInfo.FinalPos3D;
FunScore=SLMPosInfo.FinalFunScore;
Group=SLMPosInfo.Group;
confSet=SLMPosInfo.confSetFinal;
Zdepth=confSet.scan_Z+confSet.ETL

for iGroup=1:length(Group)
    Pos3DGroup{iGroup}=Pos3DAll(Group(iGroup).Indices,:);
end

PSTHparam.PreSLMCal=15;        %<----------------------------------------------------------------------------------Edit,Frame # before SLM to calculate baseline map
PSTHparam.PostSLMCal=10;        %<----------------------------------------------------------------------------------Edit,Frame # after SLM to calculate responsive map
PSTHparam.YLim=[-50 600];       % %<-----------------------------------------For Suite2p based ROI signal only method
PSTHparam.pTh=0.05;             % %<-----------------------------------------For Suite2p based ROI signal only method
PSTHparam.TestMethod='ranksum'; % %<-----------------------------------------For Suite2p based ROI signal only method
PSTHparam.FrameStep=3;          %%<----------------------------------------- Edit, time window size for dyanmic changes of PSTH.
PSTHparam.MPFrameJump=2;


idRanges=[37;58]

idRanges=[59;65]

DataList=dir([DataFolder,'ExpInfo*.mat']);

iFile=1;
A=load([DataFolder DataList(iFile).name])
TBLmat=[A.FileGenerateInfo.FileID, A.FileGenerateInfo.motionMed];

[MatFile, MatExp] = ExtractExp_FromMat(DataFolder);
MatFile.FileType=zeros(size(MatFile,1),1);
MatFile.FileType(MatFile.FileID>=55&MatFile.FileID<=64)=1;
MatFile.FileType(MatFile.FileID>=65)=2;

idRanges=repmat(MatFile.FileID(MatFile.FileType>0&MatFile.motionMed<10)',2,1);
[PSTHall, OutTBLAll] = getSLMGroup_BinNonMat(DataFolder, confSet, PSTHparam, Pos3DGroup, idRanges);
OutTBLorigin=OutTBLAll;

OutTBLorigin = join(OutTBLorigin, MatFile, 'Keys', 'FileID');

load('C:\Users\zhangl33\Projects\Project-SLMonlineControl\config\PreGenerateTseriesMultiZ\SpontBeh5T_Z11Frame550.mat','TSeriesBrukerTBL');
TSeriesBrukerTBL1=TSeriesBrukerTBL;
load('C:\Users\zhangl33\Projects\Project-SLMonlineControl\config\PreGenerateTseriesMultiZ\Anesthesia5T_Z11Frame550.mat','TSeriesBrukerTBL');
TSeriesBrukerTBL2=TSeriesBrukerTBL;
clear TSeriesBrukerTBL

TSeriesBrukerTBL=[TSeriesBrukerTBL1 TSeriesBrukerTBL2];
OutTBLorigin = MatchOutTBLAll_TSeriesBruker(OutTBLorigin, TSeriesBrukerTBL);

OutTBLAll=OutTBLorigin; 
OutTBLAll.AwakeState(OutTBLAll.TSeriesInd<=5)=1;
OutTBLAll.AwakeState(OutTBLAll.TSeriesInd>=6)=2;



MovieTable=OutTBLAll(OutTBLAll.FileID==61,:)

FrameLim=[280;360];
InvalidFrame=MovieTable.markCycle(MovieTable.markCycle>=FrameLim(1)&MovieTable.markCycle<=FrameLim(2));
GroupSeq=MovieTable.Group(MovieTable.markCycle>=FrameLim(1)&MovieTable.markCycle<=FrameLim(2))

InvalidFrame=union(InvalidFrame,InvalidFrame+1);

FrameNeed=setdiff(FrameLim(1):FrameLim(2),InvalidFrame);
Break=find(diff(FrameNeed)>1);
Break=[1 Break length(FrameNeed)];
StateSE=[Break(1:end-1);Break(2:end)];
FrameState=zeros(1,length(FrameNeed));
for iS=1:size(StateSE,2)
    FrameState(StateSE(1,iS):StateSE(2,iS))=iS; 
end


ShowMPframe=zeros(size(FrameState));
ShowMPStepAdvance=15;
ShowMPStepPost=9;
ShowMPframePost=zeros(size(FrameState));

for iS=1:size(StateSE,2)-1
    ShowMPframe(StateSE(2,iS)-ShowMPStepAdvance:StateSE(2,iS))=GroupSeq(iS);
    ShowMPframePost(StateSE(2,iS)+1:StateSE(2,iS)+ShowMPStepPost)=GroupSeq(iS);
end

% figure;
% plot(ShowMPframe);
% hold on;
% plot(ShowMPframePost)


BinFile=dir([DataFolder,'Tseries*' num2str(MovieTable.FileID(1)) '.bin']);
BinFile=[DataFolder BinFile.name];
[Data3Plane,ValidFrame] = Suite2pSingleChBin2Frame(BinFile, confSet.SLM_Pixels_Y, confSet.SLM_Pixels_X, length(confSet.ETL), FrameNeed);


nPlane=length(confSet.ETL)
for iPlane=1:nPlane
    DataPlane{iPlane}=SmoothDecDim3(squeeze(Data3Plane(:,:,:,iPlane)),1);
end

windowSize=3;
% ImgSmooth = smoothImageSequence(DataPlane{1}, FrameState, windowSize);
% ImgSmooth = ImageSequenceSlidingWin(DataPlane{1}, windowSize);
clear Data3PlaneSmooth
for iPlane=1:nPlane
    Data3PlaneSmooth(:,:,:,iPlane)= ImageSequenceSlidingWin(DataPlane{iPlane}, windowSize);
end

NormDim=4;
prcTh=[2 99.5];

Data3PlaneSmooth=AmpNormalizeDim(Data3PlaneSmooth,NormDim,prcTh);


figure;
Zdepth=unique(Pos3DAll(:,3))
ImgClim=[0.05 1];
colorCell=[0 1 0];
        PlotParam.RowPlot=1;
        PlotParam.RowColNum=1;
        PlotParam.RowColID=1;
        PlotParam.EdgeParam=[0.02 0.05 0.1 0.02 0.02 0.02];
        PlotParam.CellCenterWith=1;
        PlotParam.CellBoundaryWidth=0.5;

ColorGroup=[247 150 111;239 109 249;121 247 111]/255;
SaveVideo=[SaveFolder 'Video\'];
mkdir(SaveVideo)

% colorMapC=colormap(gray);
% colorMapC(end,:)=[1 0 0];
colorMapC=colormap(gray);

% colorMapC(end+1:end+21,:)=colorMapPN1(end-20:end,:);
close all
deleteFile(SaveVideo,'')

% figure;
for iFrame=1:size(Data3Plane,3)
% for iFrame=1:90
    figure;
    if ShowMPframe(iFrame)~=0&&ShowMPframePost(iFrame)==0
       PlotParam.CellCenterWith=2;
       H=MultiPlanes2DShow(squeeze(Data3PlaneSmooth(:,:,iFrame,:)), [], Pos3DGroup{ShowMPframe(iFrame)}, [], Zdepth, ColorGroup(ShowMPframe(iFrame),:), ImgClim,PlotParam);
       title(H,['SLM target Group ' num2str(ShowMPframe(iFrame))],'Color',ColorGroup(ShowMPframe(iFrame),:))
    elseif ShowMPframe(iFrame)==0&&ShowMPframePost(iFrame)~=0
       PlotParam.CellCenterWith=2.5;
       % H=MultiPlanes2DShow(squeeze(Data3PlaneSmooth(:,:,iFrame,:)), [], Pos3DGroup{ShowMPframePost(iFrame)}, [], Zdepth, ColorGroup(ShowMPframePost(iFrame),:), ImgClim,PlotParam);
       H=MultiPlanes2DShow(squeeze(Data3PlaneSmooth(:,:,iFrame,:)), [], Pos3DGroup{ShowMPframePost(iFrame)}, [], Zdepth, [1 0 0], ImgClim,PlotParam);
       title(H,'SLM applied','Color',[1 0 0]);
       set(H,'box','on','xcolor',[1 0 0],'ycolor',[1 0 0],'linewidth',3)
    else
       H=MultiPlanes2DShow(squeeze(Data3PlaneSmooth(:,:,iFrame,:)), [], [], [], Zdepth, ColorGroup(1,:), ImgClim,PlotParam);
    end
    b=colorbar;
    set(b,'position',[0.97 0.3 0.01 0.4],'ticks',ImgClim,'ticklabels',{'Low','High'})
    b.Label.String='Norm. F'

    colormap(colorMapC);
     % pause(0.1)
      papersizePX=[0 0 36 13];
      set(gcf, 'PaperUnits', 'centimeters');
      set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
      saveas(gcf,[SaveVideo num2str(iFrame)],'jpeg'); 
      close all;
end

% implay(mat2gray(ImgSmooth),6.9/2)  % Play the smoothed image sequence

fps=6.9*0.75;
createVideoFromJpeg(SaveVideo,fps)
% deleteFile(SaveVideo,'.jpg')

