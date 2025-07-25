function [GgraphOut,Res,r,p]=ResPosNetwork_ScoreNodeOrder_SinglePlots(SLMPosInfo,PowerTestAdj,PowerTestNode,ParamNet,ScoreCell)

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
t1=subplotPosLu(xMat,yMat,1,2,P);
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
bar1.Position=[t1.Position(1)+0.05,t1.Position(2)+t1.Position(4)+0.01, 0.1,0.02];
bar1.Label.String='Responses';
bar1.Ticks=union(ResponseLim,0);

t2=subplotPosLu(xMat,yMat,1,1,P);
imagesc([ScoreCell(rankI)]);
colormap(t2,ScoreMap);
xticklabels({'Speed','WiskStim'})
set(gca,'clim',ScoreLim,'xtick',[],'ytick',[]);
hold on;
for iFun=1:length(SLMPosInfo.Group)
    plot(0.5,TargetCellSortedP{iFun},'o','Color',GroupColor(iFun,:),'MarkerFaceColor',GroupColor(iFun,:));
end




% ClimLink=[-0.1 0.1];
% ClimCorr=[-0.4 0.4];

radius=5;
b=PowerTestAdj(rankI,rankI);
b(abs(b)<0.01)=NaN;
a=abs(b);


PosRate=sum(AdjM.*(AdjM>0))';
NegRate=sum(AdjM.*(AdjM<0))';

AllRate=sum(AdjM)';

t3=subplotPosLu(xMat,yMat,1,3,P)
[GgraphOut.G,GgraphOut.p]=MarkovState_HeatStrPlot(a,PowerTestNode(rankI)~=0,b,PosRate(rankI),NodeColorSorted,ResponseMap,ResponseLim);
GgraphOut.p.NodeColor=NodeColorSorted;
theta=ClockWiseGraph(GgraphOut.G,GgraphOut.p,radius);
% p.MarkerSize=p.MarkerSize;

CircosBand(1).rband=[radius+0.5;radius+2];
CircosBand(1).Values=PosRate(rankI);
CircosBand(1).Clim=ResponseLim;
CircosBand(1).theta=theta;
CircosBand(1).Colormap=ResponseMap;
CircosBand(1).arc_gap=abs(theta(2)-theta(1))/10;

CircosBand(end+1).rband=CircosBand(end).rband+1.2;
CircosBand(end).Values=NegRate(rankI);
CircosBand(end).Clim=ResponseLim;
CircosBand(end).theta=theta;
CircosBand(end).Colormap=ResponseMap;
CircosBand(end).arc_gap=abs(theta(2)-theta(1))/10;


CircosBand(end+1).rband=CircosBand(end).rband+1.2;
CircosBand(end).Values=ScoreCell(rankI,1);
CircosBand(end).Clim=ScoreLim;
CircosBand(end).theta=theta;
CircosBand(end).Colormap=ScoreMap;
CircosBand(end).arc_gap=0;

% CircosBand(end+1)=CircosBand(end);
% CircosBand(end).rband=CircosBand(end-1).rband+1.2;
% CircosBand(end).Values=rStim(rankI,1);
% CircosBand(end).Clim=CircosBand(end-1).Clim/2;

AddCircosBand(CircosBand);
p.ArrowSize=10;

colormap(t3,ScoreMap);
set(gca,'clim',ScoreLim);
bar2=colorbar(t3);
bar2.Location='northoutside';
bar2.Position=[t3.Position(1)+0.05,t3.Position(2)+t3.Position(4)+0.01, 0.1,0.02];
bar2.Label.String=ParamNet.ScoreLabel;
bar2.Ticks=union(ScoreLim,0);



   ResParam.Color=[0.1 0.1 0.1];
   ResParam.Marker='o';
   ResParam.MarkerSize=6;
   ResParam.Rtype='spearman';
   ResParam.xLim=ScoreLim;
   ResParam.yLim=[0 ResponseLim(2)];
   ResParam.xLabel=ParamNet.ScoreLabel;
   ResParam.yLabel='Excitation';



t4=subplotPosLu(xMat,yMat,1,4,P);
if ParamNet.AllCellRes
[Res,r,p]=LuPairRegressPlot(ScoreCell(:),PosRate(:),ResParam);
else
[Res,r,p]=LuPairRegressPlot(ScoreCell(PosRate>0),PosRate(PosRate>0),ResParam);
end

set(gca,'ytick',ResParam.yLim,'xtick',union(0,ResParam.xLim))

   ResParam.Color=[0.1 0.1 0.1];
   ResParam.yLabel='Inhibition';
   ResParam.yLim=[ResponseLim(1) 0];

t5=subplotPosLu(xMat,yMat,1,5,P);
if ParamNet.AllCellRes
[Res(2),r(2),p(2)]=LuPairRegressPlot(ScoreCell(:),NegRate(:),ResParam);
else
[Res(2),r(2),p(2)]=LuPairRegressPlot(ScoreCell(NegRate<0),NegRate(NegRate<0),ResParam);
end
set(gca,'ytick',ResParam.yLim,'xtick',union(0,ResParam.xLim))


   ResParam.Color=[0.1 0.1 0.1];
   ResParam.yLabel='TotalResponse';
   ResParam.yLim=[ResponseLim(1) ResponseLim(2)];

t6=subplotPosLu(xMat,yMat,1,6,P);
if ParamNet.AllCellRes
[Res(2),r(2),p(2)]=LuPairRegressPlot(ScoreCell(:),AllRate(:),ResParam);
else
[Res(2),r(2),p(2)]=LuPairRegressPlot(ScoreCell(AllRate~=0),AllRate(AllRate~=0),ResParam);
end
set(gca,'ytick',union(0,ResParam.yLim),'xtick',union(0,ResParam.xLim))


