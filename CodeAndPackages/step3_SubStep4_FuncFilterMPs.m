%% Functional Filter of ROIs for SLM targets 
% Including top cells highly correlated associated with speed.
XTimesStdTh = 3;
MinInterVal = 20;
maxLag = 10;
Continueloop=1;
ExcludeStimInd=[1:2000 2800:4100];  
ExcludeStimInd=union(ExcludeStimInd,ExcludeStimInd+4100);


deltaFoF=F2deltaFoF(double(CaData.F),double(CaData.Fneu),double(confSet.fs));
NeuroData=AmpNormalizeRow(deltaFoF',[0 100])';
L=min(size(fSpeed,1),size(NeuroData,1))
NeuroData=NeuroData(1:L,:);
NeuroDataCell=NeuroData(:,iscell);

TopCellNRaw=TopCellN;


for iPlane=1:numPlanes
    I1=find(CaData.CellPlaneID==iPlane);
    [rSpeed(I1),pSpeed(I1)]=corr(NeuroData(ExcludeStimInd,iscell(I1)),fSpeed(ExcludeStimInd,iPlane),'type','Spearman','rows','pairwise');
end

clear rStim c
for iCell = 1:size(iscell, 1)
    iCell;
    iPlane=CaData.CellPlaneID(iCell);
    [c(:,iCell), lags] = xcorr(NeuroData(:,iscell(iCell)), fStim(:,iPlane), maxLag, 'coeff');
    PostI = find(lags >= 0);
    [~, i1] = max(abs(c(PostI,iCell)));
    rStim(iCell, 1) = c(PostI(i1),iCell);
end


[~, ~, ~, cellBoundary] = Suite2pCellIDMapFromStat(CaData.statCell, [confSet.SLM_Pixels_X confSet.SLM_Pixels_Y]);
[~,cellPlane]=ismember(Pos3D(:,3),PlaneZ);

rSpeed;
ClimRspeed=[-0.3 0.3];
ClimRstim=[-0.1 0.1];
% CellSpeedColors = valueToColor(rSpeed, ClimRspeed, slanCM('wildfire',64));
% CellStimColors = valueToColor(rStim, ClimRstim, slanCM('wildfire',64));

CellSpeedColors = valueToColor(rSpeed, ClimRspeed, parula(64));
CellStimColors = valueToColor(rStim, ClimRstim, parula(64));



ImgClim=[0 800];
        PlotParam.RowPlot=1;
        PlotParam.RowColNum=1;
        PlotParam.RowColID=1;
        PlotParam.EdgeParam=[0.06 0.1 0.1 0.06 0.06 0.06];
        PlotParam.CellCenterWith=1;
        PlotParam.CellBoundaryWidth=1;
GroupColor=[255 51 153;91 20 212;121 247 111]/255;
figure;      
H=MultiPlanes2DShow(permute(CaData.PlaneMeanImg, [2, 1, 3]), [], Pos3Dneed, [], PlaneZ, GroupColor(FunScore(:,1),:), ImgClim,PlotParam);
subplot('position',[0.95 0.3 0.01 0.3])
MultibarPlot(1:64,gray(64),0,0,1,1);
set(gca,'xlim',[0 64],'xtick',[0 32 64],'xticklabel',sort([ClimRspeed 0]),'ytick',[]);
xlabel('Speed-Corr')
camroll(90);
papersizePX=[0 0 10*numPlanes 10];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
saveas(gcf,[SavePathAllPoint 'SpeedCorr'],'fig')
saveas(gcf,[SavePathAllPoint 'SpeedCorr.png'],'png')
print(gcf,[SavePathAllPoint 'SpeedCorr.svg'], '-dsvg', '-painters')



% SaveCorrFigure(NeuroData, rSpeed, nanmean(fSpeed,2), 1:size(NeuroData,1), 'DeltaF', SavePathAllPoint ,'Speed');
close all

figure
H=MultiPlanes2DShow(permute(CaData.PlaneMeanImg, [2, 1, 3]), cellBoundary, Pos3D, [], PlaneZ, CellStimColors, ImgClim,PlotParam);
subplot('position',[0.95 0.3 0.01 0.3])
MultibarPlot(1:64,parula(64),0,0,1,1);
set(gca,'xlim',[0 64],'xtick',[0 64],'xticklabel',sort([ClimRstim]),'ytick',[]);
xlabel('Speed-Corr')
camroll(90);
papersizePX=[0 0 10*numPlanes 10];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
saveas(gcf,[SavePathAllPoint 'StimCorr'],'fig')
saveas(gcf,[SavePathAllPoint 'StimCorr.png'],'png')
print(gcf,[SavePathAllPoint 'StimCorr.svg'], '-dsvg', '-painters')

% SaveCorrFigure(NeuroData, rSpeed, nanmean(fStim,2), 1:size(NeuroData,1), 'DeltaF', SavePathAllPoint ,'Stim');
close all






figure;
subplot(1,2,1)
% hist(rSpeed,[-1 1])
DisX=[-0.4:0.02:0.4];
m=histPlotLU(rSpeed,DisX,[0.2 0.2 0.2],0.5);
xlabel('SpeedCorr');
ylabel('Cell Counts')
subplot(1,2,2)
% hist(rSpeed,[-1 1])
DisX=[0:0.01:0.3];
m=histPlotLU(rStim,DisX,[0.2 0.2 0.2],0.5);
xlabel('StimCorr')
papersizePX=[0 0 20 10];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
saveas(gcf,[SavePathAllPoint 'DisSpeedStimCorr'],'fig')
saveas(gcf,[SavePathAllPoint 'DisSpeedStimCorr.png'],'png')

close all

rSpeed=rSpeed(:);
rStim=rStim(:);

[~,rIStim]=sort(rStim,'descend');  
[~,~,rankStim]=intersect(1:numPoint,rIStim);


[~,rISpeed]=sort(rSpeed,'descend'); 
[~,~,rankSpeed]=intersect(1:numPoint,rISpeed);

%%Cell locates close to edge of the view, were not considered as SML targets.

while Continueloop
    rCenterIStim=intersect(CenterCloseI,rIStim(1:TopCellN));
    rCenterISpeed=intersect(CenterCloseI,rISpeed(1:TopCellN));

    AllFunctionI=union(rCenterIStim,rCenterISpeed);
    ColFunctionI=intersect(rCenterIStim,rCenterISpeed);
    NonFunctionI=setdiff(CenterCloseI,AllFunctionI);

    disp([TopCellN length(rCenterISpeed) length(rCenterIStim) length(ColFunctionI)])

    rCenterIStim=setdiff(rCenterIStim,ColFunctionI);
    rCenterISpeed=setdiff(rCenterISpeed,ColFunctionI);

%     disp([TopCellN length(rCenterISpeed) length(rCenterIStim) length(ColFunctionI)])

    if length(rCenterIStim)<TopCellNRaw||length(rCenterISpeed)<TopCellNRaw
        TopCellN=TopCellN+1;
        Continueloop=1;
    else
        Continueloop=0;
    end
end

%%In case too many cells, just choose the TopCellNRaw
if length(rCenterISpeed)>TopCellNRaw
   [~,temp1]=sort(rSpeed(rCenterISpeed),'descend');
   temp2=rCenterISpeed(temp1);
   rCenterISpeed=temp2(1:TopCellNRaw);
%    rSpeed(rCenterISpeed);
end
if length(rCenterIStim)>TopCellNRaw
   [~,temp1]=sort(rStim(rCenterIStim),'descend');
   temp2=rCenterIStim(temp1);
   rCenterIStim=temp2(1:TopCellNRaw);
%    rStim(rCenterIStim);
end

if length(NonFunctionI)<TopCellNRaw
   NonFunctionI=NonFunctionI;
else
   NonFunctionI=NonFunctionI(randperm(length(NonFunctionI),TopCellNRaw));
end

%% 
% plot(lags,c(:,rCenterIStim))
% figure;
% plot(NeuroDataCell(:,rCenterIStim));
% hold on;
% plot(fStim(:,1),'color',[0.4 0.4 0.4])


IncludeCellFunFilter=[rCenterISpeed(:);rCenterIStim(:);NonFunctionI(:)];
IncludeCellFunID=[repmat(1,length(rCenterISpeed),1);repmat(2,length(rCenterIStim),1);repmat(3,length(NonFunctionI),1)];
IncludeCellFunSore=[rSpeed(IncludeCellFunFilter) rStim(IncludeCellFunFilter)];

IncludeInfo=[IncludeCellFunFilter IncludeCellFunID IncludeCellFunSore];
[~,sortI]=sort(IncludeCellFunFilter);
IncludeInfo=IncludeInfo(sortI,:);

IncludeCellFunFilter=IncludeInfo(:,1);
FunScore=IncludeInfo(:,2:end);
SavePathStimSpeed=[SavePath 'Top' num2str(TopCellN) 'SpeedStimEdgeExc\']
mkdir(SavePathStimSpeed)


IncludePathOri=[SavePathStimSpeed '\AllIncludedOrigin\'];
IncludePath=[SavePathStimSpeed '\AllIncluded\'];
mkdir(IncludePathOri)
mkdir(IncludePath)




TempGroup={rCenterISpeed(:) rCenterIStim(:) NonFunctionI(:)}

GroupName={'L cell','S cell','Non-LS cell'}
for iGroup=1:length(TempGroup)
figure;      
H=MultiPlanes2DShow(permute(CaData.PlaneMeanImg, [2, 1, 3]), cellBoundary(TempGroup{iGroup}), Pos3D(TempGroup{iGroup},:), [], PlaneZ, CellSpeedColors(TempGroup{iGroup},:), ImgClim,PlotParam);
subplot('position',[0.95 0.3 0.01 0.3])
MultibarPlot(1:64,jet(64),0,0,1,1);
set(gca,'xlim',[0 64],'xtick',[0 32 64],'xticklabel',sort([ClimRspeed 0]),'ytick',[]);
xlabel('Speed-Corr')
camroll(90);
papersizePX=[0 0 10*numPlanes 10];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
saveas(gcf,[IncludePathOri GroupName{iGroup} 'SpeedCorr'],'fig')
saveas(gcf,[IncludePathOri GroupName{iGroup} 'SpeedCorr.png'],'png')
end
close all
for iGroup=1:length(TempGroup)
figure;      
H=MultiPlanes2DShow(permute(CaData.PlaneMeanImg, [2, 1, 3]), cellBoundary(TempGroup{iGroup}), Pos3D(TempGroup{iGroup},:), [], PlaneZ, CellStimColors(TempGroup{iGroup},:), ImgClim,PlotParam);
subplot('position',[0.95 0.3 0.01 0.3])
MultibarPlot(1:64,jet(64),0,0,1,1);
set(gca,'xlim',[0 64],'xtick',[0 64],'xticklabel',sort([ClimRstim]),'ytick',[]);
xlabel('Speed-Corr')
camroll(90);
papersizePX=[0 0 10*numPlanes 10];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
saveas(gcf,[IncludePathOri GroupName{iGroup} 'StimCorr'],'fig')
saveas(gcf,[IncludePathOri GroupName{iGroup} 'StimCorr.png'],'png')
end
close all


XYZtoMarkPoint(IncludePathOri,Pos3D,IncludeCellFunFilter,yaml,confSet,CaData.statCell,FunScore);






