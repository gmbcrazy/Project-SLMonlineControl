function ProcessFOVGroupComparison(rSpeed, rStim, TargetCellList, TargetCellListFunGroup, SLMPosInfo, Nlabel, ResultFolder)
% Group comparison and boxplot for FOV

nGroup = length(SLMPosInfo.Group);
GroupLabel = {'L.','S.','N.'}; % Adjust as needed
GroupColor = [247 150 111;239 109 249;121 247 111]/255;

for iData = 1:length(Nlabel)
    Data1 = cell(1,2*nGroup);
    Data2 = cell(1,2*nGroup);
    for iFun = 1:nGroup
        Data1{iFun} = rSpeed(TargetCellList(TargetCellListFunGroup == iFun),1,iData);
        Data1{iFun + nGroup} = rSpeed(TargetCellList(TargetCellListFunGroup == iFun),2,iData);
        Data2{iFun} = rStim(TargetCellList(TargetCellListFunGroup == iFun),1,iData);
        Data2{iFun + nGroup} = rStim(TargetCellList(TargetCellListFunGroup == iFun),2,iData);
    end
    x1 = [1:nGroup, [1:nGroup] + nGroup];
    x2 = [1:nGroup, [1:nGroup] + nGroup];

    GroupPair.CorrName = 'fdr';
    GroupPair.Q = 0.1;
    GroupPair.SignY = 1;
    GroupPair.Plot = 1;
    GroupPair.Std = 1;
    GroupPair.SamplePlot = 1;
    GroupPair.SamplePairedPlot = 1;
    GroupPair.LimY = [0 GroupPair.SignY*1.2];
    GroupPair.Marker = {'o'};
    GroupPair.GroupName = repmat(GroupLabel,1,2);
   p1=[1 1 2;2 3 3];p1=[p1 p1+nGroup];
   p2=[1:nGroup];p2=[p2;p2+nGroup];
   GroupPair.Pair=[p1 p2];
    GroupPair.GroupName=repmat(GroupLabel,1,2);
   GroupPair.SignY=0.5;

    % P.xLeft=0.06;        %%%%%%Left Margin
    P.xRight=0.02;       %%%%%%Right Margin
    P.yTop=0.02;         %%%%%%Top Margin
    P.yBottom=0.1;      %%%%%%Bottom Margin
    % P.xInt=0.02;         %%%%%%Width-interval between subplots
    P.yInt=0.01;         %%%%%%Height-interval between subplots
    P.xInt = 0.1;
    P.xLeft = 0.1;
    
    figure;
    subplotLU(1,2,1,1,P);
    ErrorBarPlotLU(x1, Data1, [], repmat(GroupColor,2,1),2,1,[ResultFolder Nlabel{iData} 'SpeedGroup.txt'],GroupPair,repmat(1:nGroup,1,2));
    LuLegend([9 9 9;0.5 0.4 0.3;10 10 10;0.5 0.4 0.3],0,GroupLabel,GroupColor,8);
    set(gca,'xlim',[0 12],'ylim',[-0.2 0.6],'ytick',-0.2:0.2:0.6,'xtick',[2 6],'xticklabel',{'Spon 1','Spon 2'});
    ylabel('Speed Corr')
    
    subplotLU(1,2,1,2,P);
    ErrorBarPlotLU(x2, Data2, [], repmat(GroupColor,2,1),2,1,[ResultFolder Nlabel{iData} 'StimGroup.txt'],GroupPair,repmat(1:nGroup,1,2));
    set(gca,'xlim',[0 12],'ylim',[-0.05 0.3],'ytick',-0.05:0.05:0.3,'xtick',[2 6],'xticklabel',{'Spon 1','Spon 2'});
    ylabel('Stim Corr')
    
    papersizePX=[0 0 22 10];
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
    saveas(gcf,[ResultFolder 'BehCorrCompare' Nlabel{iData}],'png');
    close all
end

end
