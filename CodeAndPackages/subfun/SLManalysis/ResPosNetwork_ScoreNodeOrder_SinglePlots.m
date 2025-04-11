function [GgraphOut,Res,r,p]=ResPosNetwork_ScoreNodeOrder_SinglePlots(SLMPosInfo,PowerTestAdj,PowerTestNode,ParamNet,ScoreCell,ResultFolder)

[~,rankI]=sort(ScoreCell,'descend');

GroupColor=ParamNet.GroupColor;
TargetCellList=ParamNet.TargetCellList;
TargetCellListFunGroup=ParamNet.TargetCellListFunGroup;
iscell=ParamNet.iscell;
SuccTarget=ParamNet.SuccTarget;
ScoreLim=ParamNet.ScoreLim;
ResponseLim=ParamNet.ResponseLim;
ScoreMap=ParamNet.ScoreMap;
ResponseMap=ParamNet.ResponseMap;
NodeColor=ParamNet.NodeColor;
PowerTargetI=ParamNet.PowerTargetI;
statCellRes=ParamNet.statCellRes;
crit_pAll=ParamNet.crit_pAll;
SuccAmp=ParamNet.SuccAmp;
xMat=ParamNet.xMat;
yMat=ParamNet.yMat;

   P.xLeft=0.03;
   P.xRight=0.02;
   P.yTop=0.15;
   P.yBottom=0.1;



TargetCellBand=zeros(length(iscell),1);
TargetCellBand(TargetCellList)=TargetCellListFunGroup;


% TargetCellListFunGroup
TargetCellSortedP={};
for iFun=1:length(SLMPosInfo.Group)
    TargetCellSortedP{iFun}=find(TargetCellBand(rankI)==iFun);
    NodeColor(TargetCellBand==iFun,:)=repmat(GroupColor(iFun,:),sum(TargetCellBand==iFun),1);
    
end
NodeColorSorted=NodeColor(rankI,:);
% TargetCellSortedP=find(TargetCellBand(rankI)>=1);



AdjM=PowerTestAdj;
AdjM(isnan(PowerTestAdj))=0;

figure;
t1=subplot('position',[0.08 0.08 0.8 0.8]);
imagesc(AdjM(rankI,rankI));
hold on;
for iFun=1:length(SLMPosInfo.Group)
    plot(0,TargetCellSortedP{iFun},'o','Color',GroupColor(iFun,:),'MarkerFaceColor',GroupColor(iFun,:))
end
colormap(t1,ResponseMap)
set(gca,'clim',ResponseLim,'xtick',[],'ytick',[],'xlim',[-0.5 size(AdjM,2)+1])
% axis off
xlabel('Cells')
ylabel('Cells')
bar1=colorbar(t1)
bar1.Location='northoutside';
bar1.Position=[0.3 0.95, 0.4,0.01];
bar1.Label.String='Responses';
papersizePX=[0 0 8 10];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
saveas(gcf,[ResultFolder ParamNet.ScoreLabel 'SubP1Res'],'tif');
saveas(gcf,[ResultFolder ParamNet.ScoreLabel 'SubP1Res.eps'],'epsc');
saveas(gcf,[ResultFolder ParamNet.ScoreLabel 'SubP1Res'],'svg');




figure;
t2=subplot('position',[0.08 0.08 0.8 0.8]);
imagesc([ScoreCell(rankI)]);
colormap(t2,ScoreMap);
xticklabels({'Speed','WiskStim'})
set(gca,'clim',ScoreLim,'xtick',[],'ytick',[]);
hold on;
for iFun=1:length(SLMPosInfo.Group)
    plot(0.5,TargetCellSortedP{iFun},'o','Color',GroupColor(iFun,:),'MarkerFaceColor',GroupColor(iFun,:));
end
papersizePX=[0 0 2 10];
set(gcf, 'PaperUnits', 'centimeters');

set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
saveas(gcf,[ResultFolder ParamNet.ScoreLabel 'SubP2SortedScore'],'tif');
saveas(gcf,[ResultFolder ParamNet.ScoreLabel 'SubP2SortedScore.eps'],'epsc');
saveas(gcf,[ResultFolder ParamNet.ScoreLabel 'SubP2SortedScore'],'svg');




% ClimLink=[-0.1 0.1];
% ClimCorr=[-0.4 0.4];

radius=5;
b=PowerTestAdj(rankI,rankI);
b(abs(b)<0.01)=NaN;
a=abs(b);


PosRate=sum(AdjM.*(AdjM>0));
NegRate=sum(AdjM.*(AdjM<0));

figure;
t3=subplot('position',[0.08 0.08 0.8 0.8]);
[GgraphOut.G,GgraphOut.p]=MarkovState_HeatStrPlot(a,PowerTestNode(rankI)~=0,b,PosRate(rankI),NodeColorSorted,ResponseMap,ResponseLim);
GgraphOut.p.NodeColor=NodeColorSorted;
theta=ClockWiseGraph(GgraphOut.G,GgraphOut.p,radius);
% p.MarkerSize=p.MarkerSize;

CircosBand(1).rband=[radius+0.5;radius+2];
CircosBand(1).Values=PosRate(rankI);
CircosBand(1).Clim=ResponseLim;
CircosBand(1).theta=theta;
CircosBand(1).Colormap=ResponseMap;
CircosBand(1).arc_gap=abs(theta(2)-theta(1))/5;

CircosBand(end+1).rband=CircosBand(end).rband+1.2;
CircosBand(end).Values=NegRate(rankI);
CircosBand(end).Clim=ResponseLim;
CircosBand(end).theta=theta;
CircosBand(end).Colormap=ResponseMap;
CircosBand(end).arc_gap=abs(theta(2)-theta(1))/5;


CircosBand(end+1).rband=CircosBand(end).rband+1.2;
CircosBand(end).Values=ScoreCell(rankI,1);
CircosBand(end).Clim=ScoreLim;
CircosBand(end).theta=theta;
CircosBand(end).Colormap=ScoreMap;
CircosBand(end).arc_gap=abs(theta(2)-theta(1))/5;

% CircosBand(end+1)=CircosBand(end);
% CircosBand(end).rband=CircosBand(end-1).rband+1.2;
% CircosBand(end).Values=rStim(rankI,1);
% CircosBand(end).Clim=CircosBand(end-1).Clim/2;

% AddCircosBand(CircosBand);
% p.ArrowSize=10;

colormap(t3,ScoreMap);
set(gca,'clim',ScoreLim);
% bar2=colorbar(t3);
% bar2.Location='northoutside';
% bar2.Position=[0.3 0.95, 0.4,0.01];
% bar2.Label.String=ParamNet.ScoreLabel;
set(gca,'xlim',radius*[-2 2],'ylim',radius*[-2 2],'xtick',[],'ytick',[]);
papersizePX=[0 0 20 20];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
% saveas(gcf,[ResultFolder ParamNet.ScoreLabel 'SubP3Circos'],'tif');
% saveas(gcf,[ResultFolder ParamNet.ScoreLabel 'SubP3Circos.eps'],'epsc');
% saveas(gcf,[ResultFolder ParamNet.ScoreLabel 'SubP3Circos'],'svg');
% 


% figure;
% t3=subplot('position',[0.08 0.08 0.8 0.8]);
AddCircosBand(CircosBand);
% p.ArrowSize=10;
% colormap(t3,ScoreMap);
% drawnow

set(gca,'clim',ScoreLim);
% bar2=colorbar(t3);
% bar2.Location='northoutside';
% bar2.Position=[0.3 0.95, 0.4,0.01];
% bar2.Label.String=ParamNet.ScoreLabel;
set(gca,'xlim',radius*[-2 2],'ylim',radius*[-2 2],'xtick',[],'ytick',[]);
papersizePX=[0 0 20 20];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
saveas(gcf,[ResultFolder ParamNet.ScoreLabel 'SubP3Circos2'],'tif');
saveas(gcf,[ResultFolder ParamNet.ScoreLabel 'SubP3Circos2.eps'],'epsc');
% saveas(gcf,[ResultFolder ParamNet.ScoreLabel 'SubP3Circos2'],'svg');
% ax=gca;
% exportgraphics(gcf, [ResultFolder ParamNet.ScoreLabel 'SubP3Circos2.pdf'], 'ContentType', 'vector');
print(gcf, [ResultFolder ParamNet.ScoreLabel 'SubP3Circos2.svg'], '-dsvg', '-painters');


   ResParam.Color=[0.1 0.1 0.1];
   ResParam.Marker='o';
   ResParam.MarkerSize=6;
   ResParam.Rtype='pearson';
   ResParam.xLim=ScoreLim;
   ResParam.yLim=[0 ResponseLim(2)];
   ResParam.xLabel=ParamNet.ScoreLabel;
   ResParam.yLabel='Excitation';



figure;
t4=subplot('position',[0.08 0.08 0.8 0.8]);
[Res,r,p]=LuPairRegressPlot(ScoreCell(:),PosRate(:),ResParam);
set(gca,'ytick',ResParam.yLim,'xtick',union(0,ResParam.xLim))
papersizePX=[0 0 10 10];
set(gcf, 'PaperUnits', 'centimeters');

saveas(gcf,[ResultFolder ParamNet.ScoreLabel 'SubP4Exc'],'tif');
saveas(gcf,[ResultFolder ParamNet.ScoreLabel 'SubP4Exc.eps'],'epsc');
saveas(gcf,[ResultFolder ParamNet.ScoreLabel 'SubP4Exc'],'svg');



   ResParam.Color=[0.1 0.1 0.1];
   ResParam.yLabel='Inhibition';
   ResParam.yLim=[ResponseLim(1) 0];

figure;
t5=subplot('position',[0.08 0.08 0.8 0.8]);
[Res(2),r(2),p(2)]=LuPairRegressPlot(ScoreCell(:),NegRate(:),ResParam);
papersizePX=[0 0 10 10];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
set(gca,'ytick',ResParam.yLim,'xtick',union(0,ResParam.xLim))
saveas(gcf,[ResultFolder ParamNet.ScoreLabel 'SubP5Inh'],'tif');
saveas(gcf,[ResultFolder ParamNet.ScoreLabel 'SubP5Inh.eps'],'epsc');
saveas(gcf,[ResultFolder ParamNet.ScoreLabel 'SubP5Inh'],'svg');

close all
