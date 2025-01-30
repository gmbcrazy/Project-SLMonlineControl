
%% Load Data
clear all
% ProcessFolder='F:\LuSLMOnlineTest\04222024\SingleP\30PixelFromEdgeExc\';
% load('C:\Users\zhangl33\Projects\GenMatCode\Plotfun\Color\colorMapPN3.mat');
ConfigFolder='F:\LuSLMOnlineTest\SL0702-G7-Ai203\09192024\';
% ConfigFolder='C:\Users\zhangl33\Projects\Project-SLMonlineControl\config\';
ConfigFile='SLMsettingG7.yml';

[~,~,~,CaData,CaDataPlane,stat,yaml,confSet]=ROIToXYZ(ConfigFolder);



% load('C:\Users\User\Project-SLMonlineControl\subfun\Color\colorMapPN3.mat');
load('C:\Users\zhangl33\Projects\Project-SLMonlineControl\subfun\Color\colorMapPN3.mat');


% ProcessFolder='F:\LuSLMOnlineTest\04222024\SingleP\30PixelFromEdgeExc\';
ProcessFolder='F:\LuSLMOnlineTest\SL0702-G7-Ai203\09192024\SingleP\Top15SpeedStimEdgeExc\';


DataFolder=[ProcessFolder 'Data\'];

load([ProcessFolder 'SLMIncludedIndFromIscell.mat'])


[cellIDMap,CellPixCount,MedCenter,cellBoundary,cellIDMapSingleUnit]=Suite2pCellIDMapFromStat(CaData.statCell,[yaml.FOVsize_PX yaml.FOVsize_PY]);

PointAll=1:size(Pos3Dneed,1);
PlaneZ=confSet.ETL+confSet.scan_Z(1);

SLMtargetIDMap=cellIDMapSingleUnit(:,:,SLMIncludedIndFromIscell);

%%
SumDataFolder=[ProcessFolder '\DataSum\'];
mkdir(SumDataFolder)
close all
PointAll=1:size(Pos3Dneed,1);


LaserPower = [1.7];
RoundCheck=[1:7]
PSTHparam.Clim=[-300 300]

step4_SubStep_CheckingSLMmap;
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
colormap(colorMapPN1)
saveas(gcf,[SumDataFolder 'D8.png'],'png');



clear SLMTable;
SLMTable(:,1)=round(1:size(Pos3Dneed,1));
SLMTable(:,2)=NaN;





% ROIparam.TotalSLMIDmap=SLMtargetIDMap;  %
ROIparam.TotalSLMPos3D=Pos3Dneed;  %%such that ROIparam.PointAll=1:size(Pos3Dneed,1)
ROIparam.PointAll=PointAll;
ROIparam.PlaneZ=PlaneZ;
ROIparam.CellSize=20;  %%by um;
ROIparam.threshold_percentage=0.3;
ROIparam.thNum=15;
ROIparam.max_distance=ceil(CellSize*2/3/umPerPixel);  %% 2/3 diameter of a cell by pixel as maximal response region-SLM center distance
ROIparam.min_region_size=10;
ROIparam.PeakTh=200;
ROIparam.min_merged_region_size=20;
ROIparam.contourMethod='perim';
ROIparam.NeighbourHfWidthPixel=20;
ROIparam.umPerPixel=mean([yaml.umPerlPixelX yaml.umPerlPixelY]);
ROIparam.Colormap=colorMapPN1;                  
ROIparam.LaserPower=confSet.UncagingLaserPower;         %%

ROIparam.LaserPower=confSet.UncagingLaserPower(1:2);         %%


ROIparam.PointTest=PointAll;               %%<<<<<Edit
minTrialN=1;                               %%<<<<<Edit
FileIDrange=[1;200];                            %%<<<<<Edit
ROIparam.Clim=[-400 400];                  %%<<<<<Edit
[SLMRes,sampleN]=SLMResponseROIMap(SLMTrialMap,SLMTrialInfo,ROIparam,minTrialN,SumDataFolder,FileIDrange);

TerminalTrialN=5;
ExcludeTrialN=2;

SLMTestParam.TerminalTrialN=5;  %Trials # to define SLM responsive cells
SLMTestParam.ExcludeTrialN=2;   %Trials # to define SLM non-responsive cells
SLMTestParam.AllLaserPower=confSet.UncagingLaserPower;   %Trials # to define SLM non-responsive cells


clear FinishedI ConfirmI
% FinishedI=zeros(size(SLMRes));
% ConfirmI=FinishedI;

[UpdateXml, SLMTable, PointsTest, XMLparam, InfoListByLaser]=step3Fun_NextSLMtest(SLMRes,sampleN,ROIparam,XMLparam,SLMTestParam,SLMTable);






PointsTest=[];
UpdateXml=0;
for iLaser=1:length(ROIparam.LaserPower)
    % [~,iLaser]=ismember(ROIparam.LaserPower(iLaser),confSet.UncagingLaserPower);
    SLMList{iLaser}=find(SLMRes(:,iLaser)==1&sampleN(:,iLaser)>=TerminalTrialN);
    PriorityConfirmList{iLaser}=find(SLMRes(:,iLaser)==1&sampleN(:,iLaser)<TerminalTrialN);
    ExcludeList{iLaser}=find(SLMRes(:,iLaser)==0&sampleN(:,iLaser)>=ExcludeTrialN);
    UnkownConfirmList{iLaser}=find(SLMRes(:,iLaser)==0&sampleN(:,iLaser)<ExcludeTrialN);

    SLMtempTable=SLMTable;
    SLMtempTable(:,3)=NaN;
    SLMtempTable(SLMList{iLaser},3)=ROIparam.LaserPower(iLaser);

    SLMTable(:,2)=min(SLMtempTable(:,[2 3]),[],2);
   %%
   % if (length(PriorityConfirmList{iLaser})+PriorityConfirmList{iLaser})/
   % 
   % 
   % end

   %% 
   if ~isempty(SLMList{iLaser})
       disp([num2str(length(SLMList{iLaser})) ' Priority Points respond to laser power' num2str(ROIparam.LaserPower(iLaser)) ' : ' num2str(SLMList{iLaser}') ]);
   end
   %%
   if ~isempty(ExcludeList{iLaser})
      disp([num2str(length(ExcludeList{iLaser})) ' Points were excluded for laser ' num2str(ROIparam.LaserPower(iLaser)),' requiring higher power : ' num2str(ExcludeList{iLaser}') ]);
   end

    %% First test the MP point with response but trial # is not enough
    if ~isempty(PriorityConfirmList{iLaser})&&UpdateXml==0
       UpdateXml=1;
       PointsTest=PriorityConfirmList{iLaser};
       XMLparam.Laser=ROIparam.LaserPower(iLaser);
       disp([num2str(length(PriorityConfirmList{iLaser})) ' Priority Points need more test with laser power' num2str(ROIparam.LaserPower(iLaser)) ' : ' num2str(PriorityConfirmList{iLaser}') ]);
       % continue
       % break
    end

    %% Then, test the MP point without response while trial # is not enough
    if ~isempty(PriorityConfirmList{iLaser})&&UpdateXml==0
       UpdateXml=1;
       PointsTest= UnkownConfirmList{iLaser};
       XMLparam.Laser=ROIparam.LaserPower(iLaser);
       disp([num2str(length(UnkownConfirmList{iLaser})) ' Unkown Points need more test with laser power' num2str(ROIparam.LaserPower(iLaser)) ' : ' num2str(UnkownConfirmList{iLaser}') ]);
       % break
    end

    if isempty(PointsTest)&&(~isempty(ExcludeList{iLaser}))
       PointsTest= UnkownConfirmList{iLaser};
       
       I1=find(confSet.UncagingLaserPower>ROIparam.LaserPower(iLaser));
       if isempty(I1)&&UpdateXml==0
          XMLparam.Laser=confSet.UncagingLaserPower(min(I1));
          UpdateXml=1;
       end
    end

end





SumDataFolder='F:\LuSLMOnlineTest\SL0702-G7-Ai203\09192024\SingleP\Top15SpeedStimEdgeExc\DataSum\';
load([SumDataFolder 'SLMresponse.mat'],'PVparam','PSTHparam','SLMTrialMap','SLMTrialInfo')

PointAll=1:size(Pos3Dneed,1);
LaserPower=[1.5 1.6 1.7];
TrialNTh=3;


   P.xLeft=0.06;        %%%%%%Left Margin
   P.xRight=0.02;       %%%%%%Right Margin
   P.yTop=0.02;         %%%%%%Top Margin
   P.yBottom=0.06;      %%%%%%Bottom Margin
   P.xInt=0.02;         %%%%%%Width-interval between subplots
   P.yInt=0.01;         %%%%%%Height-interval between subplots

NeighbourRadPixel=20;
NeighbourRange=[-1 1]*NeighbourRadPixel;

PSTHparam.Clim=[-500 500]

CellSize=20;  %%by um;
threshold_percentage=0.3;
thNum=15;
max_distance=ceil(CellSize*2/3/umPerPixel);  %% 2/3 diameter of a cell by pixel as maximal response region-SLM center distance
min_region_size=10;
PeakTh=200;
min_merged_region_size=20;
TrialNTh=2;
contourMethod='perim';

clear result sampleN

figure;
for iPoint=1:length(PointAll)
    for iLaser=1:length(LaserPower)
        Point=PointAll(iPoint);
        I1=find((SLMTrialInfo(:,2)==Point&ismember(SLMTrialInfo(:,3),LaserPower(iLaser)))==1);
        I2=find(SLMTrialInfo(:,1)<40);
        I1=intersect(I1,I2);

        if length(I1)>=TrialNTh          
           subplotLU(length(PointAll),length(LaserPower),iPoint,iLaser,P);
           PlaneI=find(abs(Pos3Dneed(Point,3)-PlaneZ)<1);
           PSTHtemp=squeeze(mean(SLMTrialMap(:,:,:,I1),4));     

           % [roiNeighbour,ROIboundary]=ROI2neighbour(squeeze(PSTHtemp(:,:,PlaneI))',squeeze(cellIDMapSingleUnit(:,:,SLMIncludedIndFromIscell(iPoint))),1,NeighbourRadPixel);

           tempMean=squeeze(CaData.PlaneMeanImg(:,:,PlaneI));
           tempMean1=AmpNormalize(tempMean,[0 100]);  

           heatMap=squeeze(PSTHtemp(:,:,PlaneI))';
           % % figure;
           % % subplot(1,2,1);
           % tempMean2=AmpNormalize(tempMean1,[0 90]);  

        % [roiImg,roiBImg]=ROI2neighbour(tempMean1,squeeze(cellIDMapSingleUnit(:,:,SLMIncludedIndFromIscell(iPoint))),1,NeighbourRadPixel);
        [roiHeat,roiBHeat]=ROI2neighbour(heatMap,squeeze(cellIDMapSingleUnit(:,:,SLMIncludedIndFromIscell(iPoint))),1,NeighbourRadPixel);

        [result(iPoint,iLaser), mergedRegion, contourPixels] = check_high_value_center(roiHeat, threshold_percentage, max_distance, min_region_size, PeakTh, min_merged_region_size, contourMethod);
           

        sampleN(iPoint,iLaser)=length(I1);

        imagesc(roiHeat);hold on;
        set(gca,'clim',[-400 400],'xtick',[],'ytick',[])
        colormap(colorMapPN1)
        % for iField=1:length(regions)
        if ~isempty(mergedRegion)
            plot(contourPixels(:,2),contourPixels(:,1),'g.');
        end
         text(NeighbourRadPixel,0,['P' num2str(Point) ', ' num2str(length(I1))],'color',[0 1 0]) 



        end


    end
end

papersizePX=[0 0 length(LaserPower)*4 length(PointAll)*4]
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
saveas(gcf,[SumDataFolder 'SummaryROI2.png'],'png');





figure;
PSTHparam.Clim=[-400 400]
for iPoint=1:length(PointAll)
    for iLaser=1:length(LaserPower)

        Point=PointAll(iPoint);
%         subplot(ShowRow,ShowCol,iPoint)
        I1=find((SLMTrialInfo(:,2)==Point&ismember(SLMTrialInfo(:,3),LaserPower(iLaser)))==1);
        if length(I1)>=TrialNTh          
           subplotLU(length(PointAll),length(LaserPower),iPoint,iLaser,P);
           PlaneI=find(abs(Pos3Dneed(Point,3)-PlaneZ)<1);
           PSTHtemp=squeeze(mean(SLMTrialMap(:,:,:,I1),4));     

           MultPlaneIs2DShow1Plane(PSTHtemp, [], Pos3Dneed(Point,:), [], PlaneZ, PlaneI, [0 1 0], PSTHparam.Clim);
           % text(confSet.SLM_Pixels_X/2,0, ['P' num2str(Point) ', n = ' num2str(length(I1))],'color',[0 1 0]);
           text(Pos3Dneed(Point,1),Pos3Dneed(Point,2)+NeighbourRadPixel/2, ['P' num2str(Point) ', n = ' num2str(length(I1))],'color',[0 1 0]);
           
           set(gca,'xtick',[],'ytick',[],'xlim',Pos3Dneed(Point,1)+NeighbourRange,'ylim',Pos3Dneed(Point,2)+NeighbourRange);

        end


    end
end
colormap(colorMapPN1)
papersizePX=[0 0 length(LaserPower)*8 length(PointAll)*8]
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
saveas(gcf,[SumDataFolder 'Summary2.png'],'png');

Point=5;
iLaser=1;
Xstep=50;
I1=find((SLMTrialInfo(:,2)==Point&ismember(SLMTrialInfo(:,3),LaserPower(iLaser)))==1);
PlaneI=find(abs(Pos3Dneed(Point,3)-PlaneZ)<1);
PSTHtrial=squeeze(SLMTrialMap(:,:,PlaneI,I1));   
PSTHtrial=permute(PSTHtrial,[2 1 3]);
% MultiMatrixImg2DSlice(PSTHtrial,Xstep)
% MultiMatrix3DHeatmap(PSTHtrial)

figure;
% MultiMatrix2DMatPlot(PSTHtrial,1:512,1:512)
P.Xticklabel=[]
P.XLabel=[]
P.Yticklabel=[]
P.YLabel=[]
P.Clim=PSTHparam.Clim/2;

MultiMatrix2DPlot(PSTHtrial,1:length(I1),P)
colormap(colorMapPN1);
% set(gca,'clim',PSTHparam.Clim)
%%

% % figure;
% % for Point=1:length(PointAll)
% %     subplot(ShowRow,ShowCol,Point)
% %     I1=find((SLMTrialInfo(:,2)==Point&ismember(SLMTrialInfo(:,3),LaserPower))==1);
% %     if ~isempty(I1)
% %          PSTHtemp=squeeze(mean(SLMTrialMap(:,:,:,I1),4));
% % 
% %          MultiMatrix3DPlotZ(PSTHtemp,PlaneZ,0.9);
% %          caxis(PSTHparam.Clim);
% %          Radius=10;
% %          colormap(PSTHparam.ColorMap);
% %          set(gca,'xlim',[0 confSet.SLM_Pixels_Y],'ylim',[0 confSet.SLM_Pixels_X],'zlim',PlaneZ([1 end]),'ztick',PlaneZ,'View',[64 24],'zDir','reverse');
% %          % plotCellCenter3D(SinglePxyz(iPoint,:), Radius, [0 1 0],1.5);
% %          if isfield(PSTHparam,'TargetPos')
% %             plotCellCenter3D(Pos3Dneed(Point,:), Radius, [0 1 0],1.5);
% %          end
% %     end
% %     text(confSet.SLM_Pixels_X/2,confSet.SLM_Pixels_Y/2,PlaneZ(1)-5,['Point ' num2str(Point) 'Laser' num2str(LaserPower)])
% % end

% PlaneI=find(abs(Pos3Dneed(Point,3)-PlaneZ)<1)
% MultPlaneIs2DShow1Plane(PSTHtemp, [], Pos3Dneed(Point,:), Point, PlaneZ, PlaneI, [0 1 0], PSTHparam.Clim)
% 
% 
% [cellIDMap,CellPixCount,MedCenter]=Suite2pCellIDMapFromStat(CaData.statCell(SLMIncludedIndFromIscell),[confSet.SLM_Pixels_X confSet.SLM_Pixels_Y]);
% cellBoundary = CellIDMap2Boundary(cellIDMap);
% 
% BinFile='F:\LuSLMOnlineTest\MouseMei03\04252024\SingleP\20PixelFromEdgeExc\Data\TSeries-04252024-0936-096R3Laser1.5GPoint7.bin';
% XMLparam.Point=7;
% XMLparam.Laser=1.5;
% XMLparam.RoundID=3;
% PSTHparam.TargetPos=Pos3Dneed(XMLparam.Point,:);
% PSTHparam.CellStat=CaData.statCell{SLMIncludedIndFromIscell(XMLparam.Point)};
% PSTHparam.PlaneID=CaData.CellPlaneID(SLMIncludedIndFromIscell(XMLparam.Point));
% 
% median(PSTHparam.CellStat.xpix)
% median(PSTHparam.CellStat.ypix)
% 
% ROI=[PSTHparam.CellStat.xpix(:) PSTHparam.CellStat.ypix(:)];
% [PSTHtemp,ROITseries]=PSTHmapSampleCal(BinFile,PSTHparam,confSet,ROI,PSTHparam.PlaneID);
% 
% 
% 
% figure;hold on;
% for i=1:size(ROITseries,2)
%     plot(i,ROITseries(:,i),'r.')
% end
% 
% errorbar(1:size(ROITseries,2), mean(ROITseries),ste(ROITseries))
% 
% figure;
% plot(mean(ROITseries),'r.')
% 
% figure;
% imagesc(squeeze(PSTHtemp(:,:,3)))
%          caxis(PSTHparam.Clim);
%          Radius=10;
%          colormap(PSTHparam.ColorMap);
%          hold on;
% plot(median(PSTHparam.CellStat.ypix),median(PSTHparam.CellStat.xpix),'go')
% plotCellBoundary3D(cellBoundary(XMLparam.Point), [],[0 1 0],1)
% 
% 
% PSTHparam.CellStat.med
% PSTHtemp=PSTHmapSampleCal(BinFile,PSTHparam,confSet);
%          PlaneZ=confSet.ETL+confSet.scan_Z;
%          figure
%          subplot(1,2,1)
%          MultiMatrix3DPlotZ(PSTHtemp,PlaneZ,0.9);
%          caxis(PSTHparam.Clim);
%          Radius=10;
%          colormap(PSTHparam.ColorMap);
%          set(gca,'xlim',[0 confSet.SLM_Pixels_Y],'ylim',[0 confSet.SLM_Pixels_X],'zlim',PlaneZ([1 end]),'ztick',PlaneZ,'View',[64 24],'zDir','reverse');
%          % plotCellCenter3D(SinglePxyz(iPoint,:), Radius, [0 1 0],1.5);
%          if isfield(PSTHparam,'TargetPos')
%             plotCellCenter3D(PSTHparam.TargetPos, Radius, [0 1 0],1.5);
%          end
% 
% 
%          subplot(1,2,2)
%          MultiMatrix3DPlotZ(PSTHtemp,PlaneZ,0.9);
%          caxis(PSTHparam.Clim);
%          Radius=10;
%          colormap(PSTHparam.ColorMap);
%          set(gca,'xlim',[0 confSet.SLM_Pixels_Y],'ylim',[0 confSet.SLM_Pixels_X],'zlim',PlaneZ([1 end]),'ztick',PlaneZ,'View',[64 24],'zDir','reverse');
%          plotCellBoundary3D(cellBoundary(XMLparam.Point), Pos3Dneed(XMLparam.Point,3),[0 1 0],1)
% 
% 
%          papersizePX=[0 0 12 20]
%          set(gcf, 'PaperUnits', 'centimeters');
%          set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
% 
% 
% P1 1.5   P2 1.4   P3 x  P4 1.5 P5 x  P6x   P7 1.5  P8x   P9x1.55 P10 1.5  
% 
% P11 x   P12 x   P13 1.5  P14x  P15 1.5   P16 1.5  P17x P18x1.55  P19 1.5
% 
% P20 1.5  P21 1.5   P22x
% 
% for iCell=1:size(Pos3Dneed,1)
% %     disp(["Point ID: ", num2str(iCell]))
%     LaserPower(iCell)=input(['Point ID: ' num2str(iCell) ' LaserPower: ']);
% end
% 
% PointList=1:size(Pos3Dneed,1);
% ResponseI=LaserPower>0;
% PointList=PointList(ResponseI);
% LaserList=LaserPower(ResponseI);
% 
% 
% FinalList=[];
% for iCell=1:length(PointList)
%     for iRound=1:5
%         FinalList=[FinalList;PointList(iCell) LaserList(iCell) iRound];
%     end
% end
% shuffledIndices = randperm(size(FinalList,1));
% 
% FinalListRandom=FinalList(shuffledIndices,:);
% 
% 
% for iTrial=1:size(FinalListRandom)
%     XMLparam.Point=FinalListRandom(iTrial,1);
%     XMLparam.Laser=FinalListRandom(iTrial,2);
%     XMLparam.RoundID=FinalListRandom(iTrial,3);
%     PSTHparam.TargetPos=Pos3Dneed(XMLparam.Point,:);
% % PV_LinkExcuteXML(ProcessFolder,RandomDelayInterval)
% % PV_LinkExcuteXML(XMLparam,PVparam,confSet)
%     PV_LinkExcuteXML(XMLparam,PVparam,confSet,PSTHparam);
% end
% 
% Series99, iTrial=48 has problem
% 
% % xmlCSV=[ProcessFolder 'xmlListCurrent.csv'];
% % PV_LinkExcuteFolder_byRound(xmlCSV,RandomDelayInterval,MaxFrame,BreakPoint,Round,[2]);
% 
% %Point 7 Only
% 
% %%disable breakP, 2 files?
% 
% 
% %%Enable breakP, 2files.  Flush loop.
% 
% 
% %Add Saline
% %%Flush once after breakP.  Flush once.
% %File 11
% 
% 
% 
% 
% 
% %Add Saline
% %%Flush once at the break point;
% %File 13
% 
% %File 17
% 
% 
% 
% %Point 6 Only
% %%Flush once at the break point;
% %File 13
% 
% %File 17
%%
[cellIDMap,CellPixCount,MedCenter,cellBoundary]=Suite2pCellIDMapFromStat(CaData.statCell,[yaml.FOVsize_PX yaml.FOVsize_PY]);


SumDataFolder=[ProcessFolder '\DataSum\'];
mkdir(SumDataFolder)
close all
PointAll=1:size(Pos3Dneed,1);

PointAll=PointsTest;


LaserPower = [1.7];
RoundCheck=[1:7]
PSTHparam.Clim=[-300 300]

step4_SubStep_CheckingSLMmap;
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
colormap(colorMapPN1)
saveas(gcf,[SumDataFolder 'D8.png'],'png');

save([SumDataFolder 'SLMresponse.mat']);


PointAll=1:size(Pos3Dneed,1);
LaserPower=[1.5 1.6 1.7];
TrialNTh=3;


   P.xLeft=0.06;        %%%%%%Left Margin
   P.xRight=0.02;       %%%%%%Right Margin
   P.yTop=0.02;         %%%%%%Top Margin
   P.yBottom=0.06;      %%%%%%Bottom Margin
   P.xInt=0.02;         %%%%%%Width-interval between subplots
   P.yInt=0.01;         %%%%%%Height-interval between subplots

NeighbourRadPixel=50;
NeighbourRange=[-1 1]*NeighbourRadPixel;
figure;
PSTHparam.Clim=[-400 400]
for iPoint=1:length(PointAll)
    for iLaser=1:length(LaserPower)

        Point=PointAll(iPoint);
%         subplot(ShowRow,ShowCol,iPoint)
        I1=find((SLMTrialInfo(:,2)==Point&ismember(SLMTrialInfo(:,3),LaserPower(iLaser)))==1);
        if length(I1)>=TrialNTh          
           subplotLU(length(PointAll),length(LaserPower),iPoint,iLaser,P);
           PlaneI=find(abs(Pos3Dneed(Point,3)-PlaneZ)<1);
           PSTHtemp=squeeze(mean(SLMTrialMap(:,:,:,I1),4));     




           MultPlaneIs2DShow1Plane(PSTHtemp, [], Pos3Dneed(Point,:), [], PlaneZ, PlaneI, [0 1 0], PSTHparam.Clim);
           % text(confSet.SLM_Pixels_X/2,0, ['P' num2str(Point) ', n = ' num2str(length(I1))],'color',[0 1 0]);
           text(Pos3Dneed(Point,1),Pos3Dneed(Point,2)+NeighbourRadPixel/2, ['P' num2str(Point) ', n = ' num2str(length(I1))],'color',[0 1 0]);
           
           set(gca,'xtick',[],'ytick',[],'xlim',Pos3Dneed(Point,1)+NeighbourRange,'ylim',Pos3Dneed(Point,2)+NeighbourRange);

        end


    end
end
colormap(colorMapPN1)
papersizePX=[0 0 length(LaserPower)*8 length(PointAll)*8]
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
saveas(gcf,[SumDataFolder 'Summary2.png'],'png');

Point=5;
iLaser=1;
Xstep=50;
I1=find((SLMTrialInfo(:,2)==Point&ismember(SLMTrialInfo(:,3),LaserPower(iLaser)))==1);
PlaneI=find(abs(Pos3Dneed(Point,3)-PlaneZ)<1);
% PSTHtrial=squeeze(SLMTrialMap(:,:,PlaneI,I1));  
PSTHtemp=squeeze(mean(SLMTrialMap(:,:,:,I1),4));     
figure;
MultPlaneIs2DShow1Plane(PSTHtemp, [], Pos3Dneed(Point,:), [], PlaneZ, PlaneI, [0 1 0], PSTHparam.Clim);
colormap(colorMapPN1);


PSTHtrial=permute(PSTHtrial,[2 1 3]);
% MultiMatrixImg2DSlice(PSTHtrial,Xstep)
% MultiMatrix3DHeatmap(PSTHtrial)

figure;
% MultiMatrix2DMatPlot(PSTHtrial,1:512,1:512)
P.Xticklabel=[]
P.XLabel=[]
P.Yticklabel=[]
P.YLabel=[]
P.Clim=PSTHparam.Clim/2;

MultiMatrix2DPlot(PSTHtrial,1:length(I1),P)
colormap(colorMapPN1);

[roiNeighbour,ROIboundary]=ROI2neighbour(MeanFieldMap,cellIDMap,cellIDNeed,PixFromMedCenter)