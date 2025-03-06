clear all

WorkFolder='E:\LuSLMOnlineTest\SL0855-Emx1G6CII-AAV9CAMKII\03042025\';
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

PSTHparam.PreSLMCal=7;        %<----------------------------------------------------------------------------------Edit,Frame # before SLM to calculate baseline map
PSTHparam.PostSLMCal=14;        %<----------------------------------------------------------------------------------Edit,Frame # after SLM to calculate responsive map
PSTHparam.YLim=[-50 600];       % %<-----------------------------------------For Suite2p based ROI signal only method
PSTHparam.pTh=0.05;             % %<-----------------------------------------For Suite2p based ROI signal only method
PSTHparam.TestMethod='ranksum'; % %<-----------------------------------------For Suite2p based ROI signal only method
PSTHparam.FrameStep=3;          %%<----------------------------------------- Edit, time window size for dyanmic changes of PSTH.
PSTHparam.MPFrameJump=2;


idRanges=[66;66]

idRanges=[39;84]

DataList=dir([DataFolder,'ExpInfo*.mat']);

% iFile=1;
% A=load([DataFolder DataList(iFile).name])
% TBLmat=[A.FileGenerateInfo.FileID, A.FileGenerateInfo.motionMed];

[MatFile, MatExp] = ExtractExp_FromMat(DataFolder);
% MatFile.FileType=zeros(size(MatFile,1),1);
% MatFile.FileType(MatFile.FileID>=55&MatFile.FileID<=64)=1;
% MatFile.FileType(MatFile.FileID>=65)=2;
% 
% idRanges=repmat(MatFile.FileID(MatFile.FileType>0&MatFile.motionMed<10)',2,1);
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
OutTBLAll.AwakeState(OutTBLAll.TSeriesInd<=5)=1;    %%The 1st half 5 Tseries is designed for awake state.
OutTBLAll.AwakeState(OutTBLAll.TSeriesInd>=6)=2;    %%The 2nd half 5 Tseries is designed for anesia state.

PSTHparamDyn=PSTHparam;
PSTHparamDyn.PostSLMCal=32;        %<----------------------------------------------------------------------------------Edit,Frame # after SLM to calculate responsive map
PSTHparamDyn.FrameStep=4;          %%<----------------------------------------- Edit, time window size for dyanmic changes of PSTH.
[PSTHallDyn, OutTBLAllDyn] = getSLMGroup_Dyn_BinNonMat(DataFolder, confSet, PSTHparamDyn, Pos3DGroup, idRanges);
OutTBLoriginDyn=OutTBLAllDyn;
OutTBLoriginDyn = join(OutTBLoriginDyn, MatFile, 'Keys', 'FileID');
OutTBLoriginDyn = MatchOutTBLAll_TSeriesBruker(OutTBLoriginDyn, TSeriesBrukerTBL);
OutTBLAllDyn=OutTBLoriginDyn; 
OutTBLAllDyn.AwakeState(OutTBLAllDyn.TSeriesInd<=5)=1;    %%The 1st half 5 Tseries is designed for awake state.
OutTBLAllDyn.AwakeState(OutTBLAllDyn.TSeriesInd>=6)=2;    %%The 2nd half 5 Tseries is designed for anesia state.




PowerGroup=[40 420];
WiskGroup=[0 1];
WiskKey={'NoWisk','Wisk'};
FunGroup=[1 2 3];
AwakeState=[1 2];
AwakeKey={'Awake','Anes'};

colorGroup=[0 1 0;1 1 0;1 0 1];
colorGroup=[0 1 0;0 1 0;0 1 0;0 1 0];

FakeColor=[0 0 0];
% close all

   P.xLeft=0.06;        %%%%%%Left Margin
   P.xRight=0.1;       %%%%%%Right Margin
   P.yTop=0.02;         %%%%%%Top Margin
   P.yBottom=0.06;      %%%%%%Bottom Margin
   P.xInt=0.02;         %%%%%%Width-interval between subplots
   P.yInt=0.02;         %%%

  PlotParam.RowPlot=0;
  PlotParam.RowColNum=length(PSTHPlot);
  PlotParam.RowColID=1;
  PlotParam.EdgeParam=[0.06 0.1 0.06 0.06 0.01 0.01];
  PlotParam.CellCenterWith=1;
  PlotParam.CellBoundaryWidth=0.5;
  PlotParam.PlotCenter=1;

clear PSTHPlotDyn PSTHPlot PSTHShowAll PSTHShowAllDyn PSTHShowKey TargetNShow
for iState=1:length(AwakeState)
    for iWisk=1:length(WiskGroup)

        Index0=OutTBLAll.UncagingLaserPower==PowerGroup(1)&OutTBLAll.AwakeState==AwakeState(iState)&OutTBLAll.VolOut==WiskGroup(iWisk);
        PSTHZero=squeeze(nanmean(PSTHall(:,:,Index0,:),3));     
        PSTHZeroDyn=squeeze(nanmean(PSTHallDyn(:,:,Index0,:,:),3));

        clear PSTHGroupDyn PSTHGroup SampleSize;
        for iGroup=1:length(Group)
             Index=OutTBLAll.UncagingLaserPower==PowerGroup(2)&OutTBLAll.AwakeState==AwakeState(iState)&OutTBLAll.VolOut==WiskGroup(iWisk)&OutTBLAll.Group==iGroup;
             SampleSize(iGroup) = sum(Index);
             PSTHGroup{iGroup}=squeeze(nanmean(PSTHall(:,:,Index,:),3));

            for iWin=1:size(PSTHallDyn,4)
                PSTHGroupDyn{iGroup}(:,:,iWin,:)=squeeze(nanmean(PSTHallDyn(:,:,Index,iWin,:),3));
            end
        end

        PSTHPlot={PSTHZero*3};
        PSTHPlot=[PSTHGroup PSTHPlot];

        TargetPlot={Pos3DAll};
        TargetPlot=[Pos3DGroup TargetPlot];

        PSTHPlotDyn={PSTHZeroDyn*3};
        PSTHPlotDyn=[PSTHGroupDyn PSTHPlotDyn];

        PSTHShowAll{iState,iWisk}=PSTHPlot;
        PSTHShowAllDyn{iState,iWisk}=PSTHPlotDyn;



        PSTHShowKey{iState,iWisk}=[AwakeKey{iState} WiskKey{iWisk}];
        TargetNShow{iState,iWisk}=[SampleSize sum(Index0)];

        clear PSTHPlot PSTHPlotDyn;


    end


end


TargetName = {'LC', 'SC', 'NC','0Power'};
ImgClim=[-600 600]

for iState=1:length(AwakeState)
    for iWisk=1:length(WiskGroup)
        PSTHPlot=PSTHShowAll{iState,iWisk};
        TargetN=TargetNShow{iState,iWisk};
        figure;
        
        for iGroup=1:length(PSTHPlot)
         
            PlotParam.RowColID=iGroup;
        
            TempColor=repmat(colorGroup(iGroup,:),size(TargetPlot{iGroup},1),1);
            if iGroup<4
            NonCell=find(isnan(FunScore(FunScore(:,1)==iGroup,2))==1);
            TempColor(NonCell,:)=repmat(FakeColor,length(NonCell),1);
            end
            H{iGroup}=MultiPlanes2DShow(SmoothDecDim3(PSTHPlot{iGroup},1), [], TargetPlot{iGroup}, [], Zdepth, TempColor, ImgClim, PlotParam);
        
        end
            colormap(colorMapPN1)
        
            b=colorbar;
            set(b,'position',[0.93 0.3 0.02 0.5],'Ticks',[ImgClim(1) 0 ImgClim(2)])
            b.Label.String='DeltaF (PostSLM-PreSLM)'
        
        for iplane=1:length(Zdepth)
            H{1}(iplane).YLabel.String = ['Plane ' num2str(iplane)] ;
        end
        for icol = 1:length(PSTHPlot)
            % Get the handle of the axes in the first row for each column
            ax = H{icol}(1);  
            % Add the title to the corresponding axes
            title(ax, TargetName{icol}, 'FontSize', 10, 'FontWeight', 'bold');
        
             H{icol}(3).XLabel.String = ['n = ' num2str(TargetN(icol))] ;
        
        end

       papersizePX=[0 0 8*length(PSTHPlot)+2 8*length(Zdepth)+1];
      set(gcf, 'PaperUnits', 'centimeters');
      set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
      saveas(gcf,[SaveFolder PSTHShowKey{iState,iWisk} num2str(PSTHparam.PostSLMCal) 'FramePostSLM'],'png'); 
      saveas(gcf,[SaveFolder PSTHShowKey{iState,iWisk} num2str(PSTHparam.PostSLMCal) 'FramePostSLM'],'fig'); 
     % 
     close all




    end
end





%%Needs edit for Dyn window heatmap plots
  PlotParam.RowPlot=0;
  PlotParam.RowColNum=size(PSTHallDyn,4);
  PlotParam.RowColID=1;
  PlotParam.EdgeParam=[0.06 0.1 0.06 0.06 0.01 0.01];
  PlotParam.CellCenterWith=1;
  PlotParam.CellBoundaryWidth=0.5;

ImgClim=[-600 600]

for iState=1:length(AwakeState)
    for iWisk=1:length(WiskGroup)
        TargetN=TargetNShow{iState,iWisk};
        
        for iGroup=1:length(PSTHShowAllDyn{iState,iWisk})
         
            PSTHPlot=PSTHShowAllDyn{iState,iWisk}{iGroup};


            figure;
 
            for iWin=1:size(PSTHallDyn,4)
                PlotParam.RowColID=iWin;
        
               TempColor=repmat(colorGroup(iGroup,:),size(TargetPlot{iGroup},1),1);
               if iGroup<4
                  NonCell=find(isnan(FunScore(FunScore(:,1)==iGroup,2))==1);
                  TempColor(NonCell,:)=repmat(FakeColor,length(NonCell),1);
               end

               H{iWin}=MultiPlanes2DShow(SmoothDecDim3(squeeze(PSTHPlot(:,:,iWin,:)),1), [], TargetPlot{iGroup}, [], Zdepth, TempColor, ImgClim, PlotParam);
            end
                    colormap(colorMapPN1)
        
            b=colorbar;
            set(b,'position',[0.93 0.3 0.02 0.5],'Ticks',[ImgClim(1) 0 ImgClim(2)])
            b.Label.String='DeltaF (PostSLM-PreSLM)'
        
        for iplane=1:length(Zdepth)
            H{1}(iplane).YLabel.String = ['Plane ' num2str(iplane)] ;
        end
        for icol = 1:size(PSTHPlot,3)
            % Get the handle of the axes in the first row for each column
            ax = H{icol}(1);  
            % Add the title to the corresponding axes
            title(ax, [TargetName{iGroup} 'Win' num2str(icol)], 'FontSize', 10, 'FontWeight', 'bold');
        
            H{icol}(1).XLabel.String = ['n = ' num2str(TargetN(iGroup))] ;
        
        end

       papersizePX=[0 0 8*size(PSTHPlot,3)+2 8*length(Zdepth)+1];
      set(gcf, 'PaperUnits', 'centimeters');
      set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
      saveas(gcf,[SaveFolder TargetName{iGroup} PSTHShowKey{iState,iWisk} num2str(PSTHparamDyn.PostSLMCal) 'FramePostSLM'],'png'); 
      saveas(gcf,[SaveFolder TargetName{iGroup} PSTHShowKey{iState,iWisk} num2str(PSTHparamDyn.PostSLMCal) 'FramePostSLM'],'fig'); 
             close all

        end

     % % 




    end
end

%%Needs edit for Dyn window heatmap plots









% [FileGenerateInfo,fileList, fileIDs] = getExpInfoFiles_NonMat(ProcessFolder, idRanges)
Index0=OutTBLAll.UncagingLaserPower==40;
PSTHZero=squeeze(nanmean(PSTHall(:,:,Index0,:),3));

clear PSTHGroupDyn PSTHGroup

for iGroup=1:length(Group)
    Index=OutTBLAll.UncagingLaserPower>0&OutTBLAll.Group==iGroup&OutTBLAll.FileType==2
    SampleSize(iGroup) = sum(Index);
    PSTHGroup{iGroup}=squeeze(nanmean(PSTHall(:,:,Index,:),3));

    % for iWin=1:size(PSTHallDyn,4)
    %     PSTHGroupDyn{iGroup,iWin}=squeeze(nanmean(PSTHallDyn(:,:,Index,iWin,:),3));
    % end
end





PSTHPlot={PSTHZero};
TargetPlot={Pos3DAll};
PSTHPlot=[PSTHGroup PSTHPlot];
TargetPlot=[Pos3DGroup TargetPlot];
TargetN=[SampleSize sum(Index0)];
TargetName = {'Locomotion C.', 'Sensory C.', 'Non C.','Zero power'};


ImgClim=[-150;150];
colorGroup=[0 1 0;1 1 0;1 0 1];
colorGroup=[0 1 0;0 1 0;0 1 0;0 1 0];

FakeColor=[0 0 0];
% close all

   P.xLeft=0.06;        %%%%%%Left Margin
   P.xRight=0.1;       %%%%%%Right Margin
   P.yTop=0.02;         %%%%%%Top Margin
   P.yBottom=0.06;      %%%%%%Bottom Margin
   P.xInt=0.02;         %%%%%%Width-interval between subplots
   P.yInt=0.02;         %%%

  PlotParam.RowPlot=0;
  PlotParam.RowColNum=length(PSTHPlot);
  PlotParam.RowColID=1;
  PlotParam.EdgeParam=[0.06 0.1 0.06 0.06 0.01 0.01];
  PlotParam.CellCenterWith=1;
  PlotParam.CellBoundaryWidth=0.5;


figure;

for iGroup=1:length(PSTHPlot)
 
    PlotParam.RowColID=iGroup;

    TempColor=repmat(colorGroup(iGroup,:),size(TargetPlot{iGroup},1),1);
    if iGroup<4
    NonCell=find(isnan(FunScore(FunScore(:,1)==iGroup,2))==1);
    TempColor(NonCell,:)=repmat(FakeColor,length(NonCell),1);
    end
    H{iGroup}=MultiPlanes2DShow(SmoothDecDim3(PSTHPlot{iGroup},1), [], TargetPlot{iGroup}, [], Zdepth, TempColor, ImgClim, PlotParam);

end
    colormap(colorMapPN1)

    b=colorbar;
    set(b,'position',[0.93 0.3 0.02 0.5],'Ticks',[ImgClim(1) 0 ImgClim(2)])
    b.Label.String='DeltaF (PostSLM-PreSLM)'

for iplane=1:length(Zdepth)
    H{1}(iplane).YLabel.String = ['Plane ' num2str(iplane)] ;
end
for icol = 1:length(PSTHPlot)
    % Get the handle of the axes in the first row for each column
    ax = H{icol}(1);  
    % Add the title to the corresponding axes
    title(ax, TargetName{icol}, 'FontSize', 10, 'FontWeight', 'bold');

     H{icol}(3).XLabel.String = ['n = ' num2str(TargetN(icol))] ;

end

papersizePX=[0 0 8*length(PSTHPlot)+2 8*length(Zdepth)+1];
      set(gcf, 'PaperUnits', 'centimeters');
      set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
      saveas(gcf,[SaveFolder num2str(PSTHparam.PostSLMCal) 'FramePostSLMresponseAna'],'png'); 
      saveas(gcf,[SaveFolder num2str(PSTHparam.PostSLMCal) 'FramePostSLMresponseAna'],'fig'); 

     close all





% Custom titles for each column

% Add titles at the top of each column

FunName = {'Locomotion C', 'Sensory C', 'Non C'}
ImgClim=[-500 500]

  PlotParam.RowPlot=0;
  PlotParam.RowColNum=size(PSTHallDyn,4);
  PlotParam.RowColID=1;
  PlotParam.EdgeParam=[0.06 0.1 0.06 0.06 0.01 0.01];
  PlotParam.CellCenterWith=1;
  PlotParam.CellBoundaryWidth=0.5;

close all
for iGroup=1:length(Group)
    figure;
 
    TempColor=repmat(colorGroup(iGroup,:),size(Pos3DGroup{iGroup},1),1);
    
    NonCell=find(isnan(FunScore(FunScore(:,1)==iGroup,2))==1);
    for iWin=1:size(PSTHGroupDyn,2)
    PlotParam.RowColID=iWin;

    TempColor(NonCell,:)=repmat(FakeColor,length(NonCell),1);
       HH{iWin}=MultiPlanes2DShow(SmoothDecDim3(PSTHGroupDyn{iGroup,iWin},1), [], Pos3DGroup{iGroup}, [], Zdepth, TempColor, ImgClim,PlotParam);
    end

    for iplane=1:length(Zdepth)
         HH{1}(iplane).YLabel.String = ['Plane ' num2str(iplane)] ;
    end
    for iWin = 1:size(PSTHGroupDyn,2)
    % Get the handle of the axes in the first row for each column
         ax = HH{iWin}(1);  
    % Add the title to the corresponding axes
         if iWin == 1
         title(ax, ['Win' num2str(iWin) ' (Frame #' num2str(PSTHparam.FrameStep) ')'], 'FontSize', 10, 'FontWeight', 'bold');
         else
         title(ax, ['Win' num2str(iWin) ], 'FontSize', 10,'FontWeight', 'normal');
             
         end
         % HH{icol}(3).XLabel.String = ['n = ' num2str(TargetN(icol))] ;
    
    end


    colormap(colorMapPN1);
    b=colorbar;
    set(b,'position',[0.93 0.3 0.02 0.5],'Ticks',[ImgClim(1) 0 ImgClim(2)]);
    b.Label.String='DeltaF (PostSLM-PreSLM)';

    papersizePX=[0 0 8*size(PSTHGroupDyn,2)+2 8*size(PSTHGroupDyn,1)+1];
      set(gcf, 'PaperUnits', 'centimeters');
      set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
      saveas(gcf,[SaveFolder FunName{iGroup} 'DynPostSLM'],'png'); 
      % saveas(gcf,[SaveFolder FunName{iGroup} 'DynPostSLM'],'fig'); 
end










for iGroup=1:length(Group)
    figure;
 
    TempColor=repmat(colorGroup(iGroup,:),size(Pos3DGroup{iGroup},1),1);
    
    NonCell=find(isnan(FunScore(FunScore(:,1)==iGroup,2))==1);

    TempColor(NonCell,:)=repmat(FakeColor,length(NonCell),1);
    MultiPlanes2DShow(SmoothDecDim3(PSTHGroup{iGroup},1), [], Pos3DGroup{iGroup}, [], Zdepth, TempColor, ImgClim,P)
    colormap(colorMapPN1)
    b=colorbar;
    set(b,'position',[0.95 0.])
end


    figure;
    MultiPlanes2DShow(SmoothDecDim3(PSTHZero,1), [], Pos3DAll, [], Zdepth, colorCell, ImgClim)
    colormap(colorMapPN1)


