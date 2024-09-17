%% Load Data
clear all
% ProcessFolder='F:\LuSLMOnlineTest\04222024\SingleP\30PixelFromEdgeExc\';
% load('C:\Users\zhangl33\Projects\GenMatCode\Plotfun\Color\colorMapPN3.mat');
ConfigFolder='C:\Users\User\Project-SLMonlineControl\config\';
% ConfigFolder='C:\Users\zhangl33\Projects\Project-SLMonlineControl\config\';

[~,~,~,CaData,CaDataPlane,stat,~,~]=ROIToXYZ(ConfigFolder);



load('C:\Users\User\Project-SLMonlineControl\subfun\Color\colorMapPN3.mat');

% ProcessFolder='F:\LuSLMOnlineTest\04222024\SingleP\30PixelFromEdgeExc\';
ProcessFolder='F:\LuSLMOnlineTest\SL0702-G7-Ai203\09162024\SingleP\Top6SpeedStimEdgeExc\';

DataFolder=[ProcessFolder 'Data\'];

load([ProcessFolder 'SLMIncludedIndFromIscell.mat'])





%% 
% RandomDelayInterval=[0 1]; %%Random delay is induced after each trial of stimulation.
% PointRepetition=1;  %%Trial Number per each xml MarkPoint stimulation.
nPlane=length(confSet.ETL);

PreMarkPointRepetition=40;
PostMarkPointRepetition=20;
frameRepetion=PreMarkPointRepetition+PostMarkPointRepetition; %%Total repepitions of Z series in T series setting.
PVparam.maxFrame=nPlane*frameRepetion;
PVparam.BreakPointFrame=PreMarkPointRepetition*nPlane;

PVparam.maxFrame=nPlane*frameRepetion;
PVparam.maxFrame=nPlane*frameRepetion;
PVparam.maxFrame=nPlane*frameRepetion;




XMLparam.ProcessFolder=ProcessFolder;
% load([XMLparam.ProcessFolder 'SLMIncludedIndFromIscell.mat']);
% PSTHparam.TargetPos=Pos3DNeed(XMLparam.Point,:);
PSTHparam.PreInd=PreMarkPointRepetition-20:PreMarkPointRepetition;
PSTHparam.PostInd=PreMarkPointRepetition+1:PreMarkPointRepetition+3;
PSTHparam.Plot=1;
PSTHparam.SmoothSD=1;
PSTHparam.ColorMap=ColorPN3;
PSTHparam.Clim=[-400 400]

% PV_LinkExcuteFolder(ProcessFolder,RandomDelayInterval,PointRepetition,MaxFrame,BreakPoint);

% Round=[1];
% ExcuteIndex=[];

% PointsTest=[4 11 17 29 30];
% 

%%
PointsTest=[1:9];  %%7 should not work
for iRound=4:5
for iPP=1:length(PointsTest)
XMLparam.Point=PointsTest(iPP);
XMLparam.Laser=1.6;
XMLparam.RoundID=iRound;
PSTHparam.TargetPos=Pos3DNeed(XMLparam.Point,:);
PSTHparam.CellStat=CaData.statCell{SLMIncludedIndFromIscell(XMLparam.Point)};

%%16 change PMT to 800
%%17 Shift FOV for match MP better. %%PMT to 780
%%19 Shift FOV for match MP better. %%PMT to 800
%%Add saline 21
%%PMT 780 25
%%PMT 800 26
%%PV power 226 229 230
%%Add Saline 56
%%Add Saline 86
%%Ch1 128 
% PV_LinkExcuteXML(ProcessFolder,RandomDelayInterval)
% PV_LinkExcuteXML(XMLparam,PVparam,confSet)
PV_LinkExcuteXML(XMLparam,PVparam,confSet,PSTHparam);
% PV_LinkExcuteXMLnoBin(XMLparam,PVparam,confSet)
end
end


% 
% 
% [cellIDMap,CellPixCount,MedCenter]=Suite2pCellIDMapFromStat(CaData.statCell(SLMIncludedIndFromIscell),[confSet.SLM_Pixels_X confSet.SLM_Pixels_Y]);
% cellBoundary = CellIDMap2Boundary(cellIDMap);
% 
% BinFile='F:\LuSLMOnlineTest\MouseMei03\04252024\SingleP\20PixelFromEdgeExc\Data\TSeries-04252024-0936-096R3Laser1.5GPoint7.bin';
% XMLparam.Point=7;
% XMLparam.Laser=1.5;
% XMLparam.RoundID=3;
% PSTHparam.TargetPos=Pos3DNeed(XMLparam.Point,:);
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
%          plotCellBoundary3D(cellBoundary(XMLparam.Point), Pos3DNeed(XMLparam.Point,3),[0 1 0],1)
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
% for iCell=1:size(Pos3DNeed,1)
% %     disp(["Point ID: ", num2str(iCell]))
%     LaserPower(iCell)=input(['Point ID: ' num2str(iCell) ' LaserPower: ']);
% end
% 
% PointList=1:size(Pos3DNeed,1);
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
%     PSTHparam.TargetPos=Pos3DNeed(XMLparam.Point,:);
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