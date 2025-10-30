function ProcessFOVCellResponsePlots(TargetCellList, TargetPointList, iscell, PVpower, PSTHparam, SuccTargetPvalue, TargetResponse, CellResponse, CellSampleN, TimBinFrame, ResultFolder, Nlabel)
% Cell response mapping and aligned plots (per target/per SLM power)

TrialThNum = 2;
responseColorMap = slanCM('seismic', 64);

P.xLeft = 0.06;
P.xRight = 0.02;
P.yTop = 0.02;
P.yBottom = 0.1;
P.xInt = 0.02;
P.yInt = 0.01;

Param.PlotType=3
Param.statisP=0;
Param.LegendShow=0;
Param.Legend=[]
% Per-cell mean traces (bar plot style)
figure;
for iCell = 1:length(TargetCellList)
    for iPower = 1:length(PVpower)
        if ~isempty(TargetResponse{iCell, iPower}) && CellSampleN(iCell, iPower) >= TrialThNum
            subplotLU(length(TargetCellList), length(PVpower), iCell, iPower, P);
            RateHist_GroupPlot(TimBinFrame+0.5, TargetResponse(iCell, iPower), [0.1 0.1 0.1],Param);
            hold on;
            text(-10,0.1,['C' num2str(TargetCellList(iCell)) 'P' num2str(TargetPointList(iCell)) 'n = ' num2str(CellSampleN(iCell, iPower))]);      
            text(-10,0.05,['p' showPvalue(SuccTargetPvalue(iCell,iPower),3)]);        
            set(gca, 'xlim', [-PSTHparam.PreSLMCal PSTHparam.PostSLMCal], 'xtick', -PSTHparam.PreSLMCal:5:PSTHparam.PostSLMCal);
            ylabel(['Target' num2str(iCell)]);
        end
    end
end
papersizePX = [0 0 length(PVpower)*4 5*length(TargetCellList) ];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', papersizePX, 'PaperSize', papersizePX(3:4));
saveas(gcf, [ResultFolder 'AllSLMTargetResponse' Nlabel], 'png');
close all

% Per-cell heatmaps (full matrix)
figure;
for iCell = 1:length(TargetCellList)
    for iPower = 1:length(PVpower)
        if ~isempty(CellResponse{iCell, iPower}) && CellSampleN(iCell, iPower) >= TrialThNum
            subplotLU(length(TargetCellList), length(PVpower), iCell, iPower, P);
            imagesc(TimBinFrame+0.5, 1:length(iscell), CellResponse{iCell, iPower});
            colormap(responseColorMap);
            set(gca, 'clim', [-0.1 0.1], 'ylim', [0 length(iscell)+1]);
            hold on;
            plot(TimBinFrame(1), TargetCellList(iCell)+0.5, 'g>');
            set(gca, 'xlim', [-PSTHparam.PreSLMCal PSTHparam.PostSLMCal], 'xtick', -PSTHparam.PreSLMCal:5:PSTHparam.PostSLMCal)
        end
    end
end
papersizePX = [0 0 length(PVpower)*3 5*length(TargetCellList) ];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', papersizePX, 'PaperSize', papersizePX(3:4));
saveas(gcf, [ResultFolder 'AllSLMTargetResponseMap' Nlabel], 'png');
close all
end
