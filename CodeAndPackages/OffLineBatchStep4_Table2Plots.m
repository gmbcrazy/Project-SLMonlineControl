ResponseMap=slanCM('seismic',64);
ScoreMap=slanCM('PiYG',64);
ScoreLim

%%Plot Cell Response vs other parameters
for iSpeedTh=1:length(PostPreDiffSpeedTh)
SaveP1=[SaveFunCon num2str(PostPreDiffSpeedTh(iSpeedTh)) '\'] ;
mkdir(SaveP1);

tbl=readtable([SaveP1 'SLMGroupResponseMinusPowerZero.csv']);
tbl=readtable([SaveP1 'SLMGroupResponse.csv']);

Indx=find(tbl.TargetCell==0&tbl.Sensory==2&tbl.Group==2);
Indy=find(tbl.TargetCell==0&tbl.Sensory==2&tbl.Group==3);

figure;
scatter(tbl.Response(Indx),tbl.Response(Indy),'Marker','.')
hold on;
plot([-0.3 0.4],[-0.3 0.4],':')

figure;
% stats=ErrorBarPlotLU([1 2],{tbl.Response(Indx),tbl.Response(Indy)},[],ProcessPar.GroupColor([1 2],:),2,1,[])


figure;
Ind=find(tbl.TargetCell==0&tbl.Sensory==1&tbl.Group~=3);
Ind=find(tbl.TargetCell==0&tbl.Sensory==1);

% Ind=find(tbl.TargetCell==0&tbl.Sensory==1);
[~,I1temp]=sort(tbl.Response(Ind),'descend');
Ind=Ind(I1temp);
% [~,I1temp]=sort(tbl.SpeedR(Ind).*tbl.TargetSpeedR(Ind),'descend');
% Ind=Ind(I1temp);
% [~,I1temp]=sort(tbl.SpeedR(Ind),'descend');
% Ind=Ind(I1temp);

theta=ClockWiseNode(length(Ind));
radius=2;
ResponseLim=[-0.2 0.2]

clear CircosBand;
CircosBand(1).rband=[radius+0.5;radius+2];
CircosBand(1).Values=tbl.Response(Ind);
CircosBand(1).Clim=ResponseLim;
CircosBand(1).theta=theta;
CircosBand(1).Colormap=ResponseMap;
CircosBand(1).arc_gap=(theta(2)-theta(1))/10;


% ScoreLim=[-0.4 0.4];
% CircosBand(end+1).rband=CircosBand(end).rband+1.6;
% % CircosBand(end).rband(2)=CircosBand(end).rband(1)+0.5;
% CircosBand(end).Values=tbl.TargetSpeedR(Ind);
% CircosBand(end).Clim=ScoreLim;
% CircosBand(end).theta=theta;
% CircosBand(end).Colormap=ScoreMap;
% CircosBand(end).arc_gap=0;

ScoreLim=[-0.6 0.6];
CircosBand(end+1).rband=CircosBand(end).rband+1.6;
CircosBand(end).Values=SmoothDec(tbl.SpeedR(Ind),10);
% CircosBand(end).Values=tbl.SpeedR(Ind);
CircosBand(end).Clim=ResponseLim;
CircosBand(end).theta=theta;
CircosBand(end).Colormap=ScoreMap;
CircosBand(end).arc_gap=(theta(2)-theta(1))/10;

ScoreLim=[-0.02 0.02];
CircosBand(end+1).rband=CircosBand(end).rband+1.6;
CircosBand(end).Values=SmoothDec(tbl.SpeedR(Ind).*tbl.TargetSpeedR(Ind),10);
% CircosBand(end).Values=tbl.SpeedR(Ind).*tbl.TargetSpeedR(Ind);
CircosBand(end).Clim=ScoreLim;
CircosBand(end).theta=theta;
CircosBand(end).Colormap=ScoreMap;
CircosBand(end).arc_gap=0;



% CircosBand(end+1)=CircosBand(end);
% CircosBand(end).rband=CircosBand(end-1).rband+1.2;
% CircosBand(end).Values=rStim(rankI,1);
% CircosBand(end).Clim=CircosBand(end-1).Clim/2;

AddCircosBand(CircosBand);





Ind=find(tbl.TargetCell==0&tbl.Sensory==2);
% Ind=find(tbl.TargetCell==0&tbl.Sensory==1);

[~,I1temp]=sort((tbl.Response(Ind)),'descend');
Ind=Ind(I1temp);
% [~,I1temp]=sort(tbl.SensoryR(Ind),'descend');
% Ind=Ind(I1temp);

theta=ClockWiseNode(length(Ind));
radius=1;
ResponseLim=[-0.2 0.2]

clear CircosBand;
CircosBand(1).rband=[radius+0.5;radius+2];
CircosBand(1).Values=tbl.Response(Ind);
CircosBand(1).Clim=ResponseLim;
CircosBand(1).theta=theta;
CircosBand(1).Colormap=ResponseMap;
CircosBand(1).arc_gap=(theta(2)-theta(1))/10;

ScoreLim=[-0.015 0.015];
CircosBand(end+1).rband(1)=CircosBand(end).rband(2)+0.1;
CircosBand(end).rband(2)=CircosBand(end).rband(1)+1;
CircosBand(end).Values=(tbl.SpeedR(Ind).*tbl.TargetSensoryR(Ind));

CircosBand(end).Clim=ResponseLim;
CircosBand(end).theta=theta;
CircosBand(end).Colormap=ScoreMap;
CircosBand(end).arc_gap=(theta(2)-theta(1))/10;

% ScoreLim=[-0.2 0.2];
% CircosBand(end+1).rband=CircosBand(end).rband+1.2;
% CircosBand(end).Values=tbl.TargetSensoryR(Ind);
% CircosBand(end).Clim=ResponseLim;
% CircosBand(end).theta=theta;
% CircosBand(end).Colormap=ScoreMap;
% CircosBand(end).arc_gap=(theta(2)-theta(1))/10;

ScoreLim=[-0.04 0.04];
CircosBand(end+1).rband=CircosBand(end).rband+1.1;
CircosBand(end).Values=(tbl.SensoryR(Ind).*tbl.TargetSpeedR(Ind));
CircosBand(end).Clim=ScoreLim;
CircosBand(end).theta=theta;
CircosBand(end).Colormap=ScoreMap;
CircosBand(end).arc_gap=0;

ScoreLim=[-0.2 0.2];
CircosBand(end+1).rband=CircosBand(end).rband+1.1;
CircosBand(end).Values=tbl.SensoryR(Ind);

CircosBand(end).Clim=ScoreLim;
CircosBand(end).theta=theta;
CircosBand(end).Colormap=ScoreMap;
CircosBand(end).arc_gap=0;

% ScoreLim=[-0.05 0.05];
% CircosBand(end+1).rband=CircosBand(end).rband+1.6;
% % CircosBand(end).rband(2)=CircosBand(end).rband(1)+0.5;
% CircosBand(end).Values=tbl.Speed(Ind);
% CircosBand(end).Clim=ScoreLim;
% CircosBand(end).theta=theta;
% CircosBand(end).Colormap=ScoreMap;
% CircosBand(end).arc_gap=0;

% CircosBand(end+1)=CircosBand(end);
% CircosBand(end).rband=CircosBand(end-1).rband+1.2;
% CircosBand(end).Values=rStim(rankI,1);
% CircosBand(end).Clim=CircosBand(end-1).Clim/2;

AddCircosBand(CircosBand);


% data1=OutResult(1).delta(:,4);
% data2=Output.rSpeed(:,1,1);
   Param.Color=ProcessPar.GroupColor;
   Param.Marker='o';
   Param.MarkerSize=12;
   Param.Rtype='pearson';
   % Param.xLim=[min(data1) max(data1)];
   % Param.yLim=[min(data2) max(data2)];
   Param.xLabel=[];
   Param.yLabel=[];
   Param.ExcludeColor=[0.5 0.5 0.5];
   Param.PlotExclude=0;
% [~,r,p]=LuPairRegressPlotGroup_ExcludeDots(data1,data2,Output.NeuroPos3DMeta(:,4),Output.GroupTargetCellAll(:,1), Param)    







% figure;
% for iWisk=1:2
%     clear TempData;
%     if iWisk==1
%        [~,rankI]=sort(Output.rSpeed(:,1,1),'descend');
%     else
%        [~,rankI]=sort(Output.rStim(:,1,1),'descend');
%     end
%     for iGroup=1:4
%         subplotLU(2,4,iWisk,iGroup)
%         imagesc(OutResult(iWisk).GroupResponse(rankI,:,iGroup));colormap(ResponseMap);set(gca,'clim',[-0.2 0.2])
%     end
% 
% end



for iWisk=1:2

    % subplot(1,2,iWisk)
    figure;
    for iGroup=1:3
    TempData{iGroup}=OutResult(iWisk).delta(NonTargetCell,iGroup)
    [~,pTFromZero(iWisk,iGroup),~,tFromZero(iWisk,iGroup)]=ttest(TempData{iGroup},0);
    [pRFromZero(iWisk,iGroup),~,RFromZero(iWisk,iGroup)]=signrank(TempData{iGroup},0);

    end
    % stats=ErrorViolinHalf([1 2 3],TempData(1:3),Param.Color,1,[SaveP1 ProcessPar.VolOutLabel{iWisk} '.txt'],GroupPair,[1 1 1]);
    stats=ErrorViolinHalf([1 2 3],TempData(1:3),Param.Color,1,[],GroupPair,[1 1 1]);
    ylabel('Cell response (Post-Pre ΔF)');
    xlabel('Stim group')
    hold on;
    plot([0 4],[0 0],'k:')
    set(gca,'ylim',[-0.05 0.1],'xlim',[0 4],'xtick',1:3,'XTickLabel',GroupPair.GroupName)
    papersizePX=[0 0 12 12];
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
    print(gcf, [SaveP1 ProcessPar.VolOutLabel{iWisk} 'DisNonTargetCellResponse.svg'], '-dsvg', '-painters');
    print(gcf, [SaveP1  ProcessPar.VolOutLabel{iWisk} 'DisNonTargetCellResponse.tif'], '-dtiffn', '-painters');





    GroupPairTarget=GroupPair;
    GroupPairTarget.SamplePairedPlot=0;
    figure;
    for iGroup=1:3
    TempData2{iGroup}=OutResult(iWisk).delta(Output.GroupTargetCellMerge{iGroup},iGroup);
    % pTemp=OutResult(iWisk).p(Output.GroupTargetCellMerge{iGroup},iGroup);
    % TempData2{iGroup}(pTemp>0.05)=[];
    % Temp
    [~,pTFromZero(iWisk,iGroup),~,tFromZero(iWisk,iGroup)]=ttest(TempData2{iGroup},0);
    [pRFromZero(iWisk,iGroup),~,RFromZero(iWisk,iGroup)]=signrank(TempData2{iGroup},0);
    end
    % stats=ErrorViolinHalf([1 2 3],TempData(1:3),Param.Color,1,[SaveP1 ProcessPar.VolOutLabel{iWisk} '.txt'],GroupPair,[1 1 1]);
    stats=ErrorViolinHalf([1 2 3],TempData2(1:3),Param.Color,0,[],GroupPair,[1 2 3]);
    ylabel('Cell response (Post-Pre ΔF)');
    xlabel('Stim group')
    hold on;
    plot([0 4],[0 0],'k:')
    set(gca,'ylim',[-0.05 1],'xlim',[0 4],'xtick',1:3,'XTickLabel',GroupPair.GroupName)
    papersizePX=[0 0 12 12];
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
    print(gcf, [SaveP1 ProcessPar.VolOutLabel{iWisk} 'DisTargetCellResponse.svg'], '-dsvg', '-painters');
    print(gcf, [SaveP1  ProcessPar.VolOutLabel{iWisk} 'DisTargetCellResponse.tif'], '-dtiffn', '-painters');




    



    % print(gcf, [SaveP1 ProcessPar.VolOutLabel{iWisk} 'DisCellResponse.svg'], '-dsvg', '-painters');
    % print(gcf, [SaveP1  ProcessPar.VolOutLabel{iWisk} 'DisCellResponse.tif'], '-dtiffn', '-painters');
    % close all
    % figure;
    % hold on;
    % for iGroup=1:4
    %     figure;
    %     histPlotLU(TempData2{iGroup},[-0.3:0.005:0.4],Param.Color(iGroup,:),0.3);
    % end





end
 




end
