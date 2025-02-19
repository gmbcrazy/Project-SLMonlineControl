function PostGPLadjusted(ProcessFolder,IncludeCellFunFilter, FunScore, yaml, confSet,Cellstat,NonTargetPath,NonTargets, IndexNonTargetTrial,CaData)

IncludePathOri=[ProcessFolder '\AllIncludedOrigin\'];
IncludePath=[ProcessFolder '\AllIncluded\'];

EditedTable=gpl2Table([IncludePath 'GPL.gpl']);
OriTable=gpl2Table([IncludePathOri 'GPL.gpl']);
% [~,~,I2]=intersect(EditedTable.Name,OriTable.Name);

[I1,I2]=ismember(EditedTable.Name,OriTable.Name);



IncludeCellFinal=IncludeCellFunFilter(I2);
if ~isempty(FunScore)
   FunScore=FunScore(I2,:);
   GroupIDAll=unique(FunScore(:,1));
   for iGroup=1:length(GroupIDAll)
       TempGroup{iGroup}=find(FunScore(:,1)==GroupIDAll(iGroup));
   end
end



SLMIncludedIndFromIscell=IncludeCellFinal;
[XYPosPixel,Z]=gplXYtoPixel(EditedTable,yaml);
Pos3Dneed=[XYPosPixel Z(:,1)];

save([IncludePath 'SLMIncludedIndFromIscell.mat'],'SLMIncludedIndFromIscell','Pos3Dneed','yaml','confSet','NonTargets','IndexNonTargetTrial','Cellstat','FunScore');

[~, ~, MedCenter, cellBoundary] = Suite2pCellIDMapFromStat(Cellstat(SLMIncludedIndFromIscell), [confSet.SLM_Pixels_X confSet.SLM_Pixels_Y]);
GroupName={'L cell','S cell','Non-LS cell'};
numPlanes=length(confSet.ETL);
ClimRspeed=[-0.2 0.2];
ClimRstim=[0 0.15];
ImgClim=[0 400];
        PlotParam.RowPlot=1;
        PlotParam.RowColNum=1;
        PlotParam.RowColID=1;
        PlotParam.EdgeParam=[0.06 0.1 0.1 0.06 0.06 0.06];
        PlotParam.CellCenterWith=1;
        PlotParam.CellBoundaryWidth=1;
CellSpeedColors = valueToColor(FunScore(:,2), ClimRspeed, jet(64));
CellStimColors = valueToColor(FunScore(:,3), ClimRstim, jet(64));
PlaneZ=confSet.ETL+confSet.scan_Z(1);
% for iGroup=1:length(TempGroup)
% figure;      
% H=MultiPlanes2DShow(permute(CaData.PlaneMeanImg, [2, 1, 3]), cellBoundary(TempGroup{iGroup}), Pos3Dneed(TempGroup{iGroup},:), [], PlaneZ, CellSpeedColors(TempGroup{iGroup},:), ImgClim,PlotParam);
% subplot('position',[0.95 0.3 0.01 0.3])
% MultibarPlot(1:64,jet(64),0,0,1,1);
% set(gca,'xlim',[0 64],'xtick',[0 32 64],'xticklabel',sort([ClimRspeed 0]),'ytick',[]);
% xlabel('Speed-Corr')
% camroll(90);
% papersizePX=[0 0 10*numPlanes 10];
% set(gcf, 'PaperUnits', 'centimeters');
% set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
% saveas(gcf,[IncludePath GroupName{iGroup} 'SpeedCorr'],'fig')
% saveas(gcf,[IncludePath GroupName{iGroup} 'SpeedCorr.png'],'png')
% end
% close all
% for iGroup=1:length(TempGroup)
% figure;      
% H=MultiPlanes2DShow(permute(CaData.PlaneMeanImg, [2, 1, 3]), cellBoundary(TempGroup{iGroup}), Pos3Dneed(TempGroup{iGroup},:), [], PlaneZ, CellStimColors(TempGroup{iGroup},:), ImgClim,PlotParam);
% % H=MultiPlanes2DShow(permute(CaData.PlaneMeanImg, [2, 1, 3]), cellBoundary(TempGroup{iGroup}), permute(Pos3Dneed(TempGroup{iGroup},:),[2 1 3]), [], PlaneZ, CellStimColors(TempGroup{iGroup},:), ImgClim,PlotParam);
% 
% subplot('position',[0.95 0.3 0.01 0.3])
% MultibarPlot(1:64,jet(64),0,0,1,1);
% set(gca,'xlim',[0 64],'xtick',[0 64],'xticklabel',sort([ClimRstim]),'ytick',[]);
% xlabel('Speed-Corr')
% camroll(90);
% papersizePX=[0 0 10*numPlanes 10];
% set(gcf, 'PaperUnits', 'centimeters');
% set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
% saveas(gcf,[IncludePath GroupName{iGroup} 'StimCorr'],'fig')
% saveas(gcf,[IncludePath GroupName{iGroup} 'StimCorr.png'],'png')
% end
% close all
% 





% NonTargetPath=

Pos3DFromGPL=[EditedTable.X EditedTable.Y EditedTable.Z];
NonTargetsTable=gpl2Table([NonTargetPath 'SelectedFromRaw.gpl']);
% NonTargetsTable=gpl2Table([NonTargetPath 'Raw.gpl']);

NonTargetsGPL=[NonTargetsTable.X NonTargetsTable.Y NonTargetsTable.Z];
EditedGPLtoMarkPoint_NT_PairGplXml(ProcessFolder, Pos3DFromGPL, SLMIncludedIndFromIscell, NonTargetsGPL, IndexNonTargetTrial, yaml, confSet,Cellstat,FunScore);
