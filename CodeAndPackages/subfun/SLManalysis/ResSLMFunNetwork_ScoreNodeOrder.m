function [GgraphOut,Res,r,p,tHandels]=ResSLMFunNetwork_ScoreNodeOrder(PSTHparam,ResponseAdj,ResponseNode,ParamNet,ScoreCell)


GroupColor=ParamNet.GroupColor;
TargetCellList=ParamNet.TargetCellList;
% TargetCellListFunGroup=ParamNet.TargetCellListFunGroup;
% SuccTarget=ParamNet.SuccTarget;
ScoreLim=ParamNet.ScoreLim;
ResponseLim=ParamNet.ResponseLim;
ScoreMap=ParamNet.ScoreMap;
ResponseMap=ParamNet.ResponseMap;
NodeColor=ParamNet.NodeColor;
% PowerTargetI=ParamNet.PowerTargetI;
statCellRes=ParamNet.statCellRes;
crit_pAll=ParamNet.crit_pAll;
% SuccAmp=ParamNet.SuccAmp;
xMat=ParamNet.xMat;
yMat=ParamNet.yMat;

TimeBinFrame=PSTHparam.TimeBinFrame;
TimeBinFrame= -PSTHparam.PreSLMCal:PSTHparam.PostSLMCal-1;


   P.xLeft=0.03;
   P.xRight=0.02;
   P.yTop=0.15;
   P.yBottom=0.1;

[~,rankI]=sort(ScoreCell,'descend');
rankI=[rankI;ParamNet.CellN+1];


TargetCellBand=zeros(ParamNet.CellN+1,1);
TargetCellBand(ParamNet.GroupTargetCell)=ParamNet.SLMGroup;
TargetCellBand(end)=ParamNet.SLMGroup;

% TargetCellListFunGroup
TargetCellSortedP={};
for iFun=1:size(GroupColor,1)
    TargetCellSortedP{iFun}=find(TargetCellBand(rankI)==iFun);
    NodeColor(TargetCellBand==iFun,:)=repmat(GroupColor(iFun,:),sum(TargetCellBand==iFun),1);
end
NodeColorSorted=NodeColor(rankI,:);
% NodeColorSorted=
% TargetCellSortedP=find(TargetCellBand(rankI)>=1);



AdjM=ResponseAdj;
AdjM(isnan(ResponseAdj))=0;
figure;
t1=subplotPosLu(xMat,yMat,1,2,P);
hold on;
if ~isempty(PSTHparam.Data)
 imagesc(TimeBinFrame+0.5,1:ParamNet.CellN,PSTHparam.Data(rankI(1:end-1),:));
 if ~isempty(TargetCellSortedP{ParamNet.SLMGroup}(1:end-1))
    plot(TimeBinFrame(1),TargetCellSortedP{ParamNet.SLMGroup}(1:end-1),'o','Color',GroupColor(ParamNet.SLMGroup,:),'MarkerFaceColor',GroupColor(ParamNet.SLMGroup,:))

 end
 set(gca,'xlim',[TimeBinFrame(1) TimeBinFrame(end)],'xtick',[-PSTHparam.PreSLMCal:5:PSTHparam.PostSLMCal],'ylim',[0 ParamNet.CellN+0.5]);
 set(gca,'clim',ParamNet.ResponseLim);colormap(ParamNet.ResponseMap);

axis ij
end
% hold on;
% for iFun=1:length(SLMPosInfo.Group)
%     plot(0,TargetCellSortedP{iFun},'>','Color',GroupColor(iFun,:),'MarkerFaceColor',GroupColor(iFun,:))
% end
colormap(t1,ResponseMap)
% set(gca,'clim',ResponseLim,'xtick',[],'ytick',[],'xlim',[-0.5 size(AdjM,2)+1])
% axis off
xlabel('Frame from SLM')
ylabel('Cells')

bar1=colorbar(t1);
bar1.Location='northoutside';
bar1.Position=[t1.Position(1)+0.05,t1.Position(2)+t1.Position(4)+0.01, 0.1,0.02];
bar1.Label.String='Responses';
bar1.Ticks=union(ResponseLim,0);

% 
t2=subplotPosLu(xMat,yMat,1,1,P);
imagesc([ScoreCell(rankI(1:end-1))]);
colormap(t2,ScoreMap);
xticklabels({'Speed','WiskStim'})
set(gca,'clim',ScoreLim,'xtick',[],'ytick',[]);
hold on;
 if ~isempty(TargetCellSortedP{ParamNet.SLMGroup}(1:end-1))
    plot(0.5,TargetCellSortedP{ParamNet.SLMGroup}(1:end-1),'o','Color',GroupColor(ParamNet.SLMGroup,:),'MarkerFaceColor',GroupColor(ParamNet.SLMGroup,:));
 end



% ClimLink=[-0.1 0.1];
% ClimCorr=[-0.4 0.4];

radius=5;
% rankI=[rankI;ParamNet.CellN+1];
b=ResponseAdj(rankI,rankI);
b(abs(b)<0.01)=NaN;
a=abs(b);

%b(TargetCellSortedP);

PosRate=sum(AdjM.*(AdjM>0));
NegRate=sum(AdjM.*(AdjM<0));

Trans=sum(AdjM.*(AdjM~=0));
a(end,TargetCellSortedP{ParamNet.SLMGroup})=NaN;
b(end,TargetCellSortedP{ParamNet.SLMGroup})=NaN;



t3=subplotPosLu(xMat,yMat,1,3,P);
[GgraphOut.G,GgraphOut.p]=MarkovState_HeatStrPlot(a,ResponseNode(rankI)~=0,b,Trans(rankI),NodeColorSorted,ResponseMap,ResponseLim);
GgraphOut.p.NodeColor=NodeColorSorted;
theta=ClockWiseGraph1(GgraphOut.G,GgraphOut.p,radius);
% p.MarkerSize=p.MarkerSize;


CircosBand(1).rband=[radius+0.5;radius+2];
CircosBand(1).Values=Trans(rankI(1:end-1));
CircosBand(1).Clim=ResponseLim;
CircosBand(1).theta=theta;
CircosBand(1).Colormap=ResponseMap;
CircosBand(1).arc_gap=(theta(2)-theta(1))/10;

% CircosBand(end+1).rband=CircosBand(end).rband+1.2;
% CircosBand(end).Values=Trans(rankI(1:end-1));
% CircosBand(end).Clim=ResponseLim;
% CircosBand(end).theta=theta;
% CircosBand(end).Colormap=ResponseMap;
% CircosBand(end).arc_gap=(theta(2)-theta(1))/10;


CircosBand(end+1).rband=CircosBand(end).rband+1.2;
CircosBand(end).Values=ScoreCell(rankI(1:end-1),1);
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
   ResParam.MarkerSize=12;
   ResParam.Rtype='spearman';
   ResParam.xLim=ScoreLim;
   ResParam.yLim=ResponseLim;
   ResParam.xLabel=ParamNet.ScoreLabel;
   ResParam.yLabel='Response';
   ResParam.ExcludeColor=GroupColor(ParamNet.SLMGroup,:);
   ResParam.PlotExclude=0;


t4=subplotPosLu(xMat,yMat,1,4,P);
[Res,r,p]=LuPairRegressPlot_ExcludeDots(ScoreCell,Trans(1:end-1)',ParamNet.OffTargetCellList,ResParam);
% set(gca,'ytick',ResParam.yLim,'xtick',union(0,ResParam.xLim))
% 
% IncludeInd=setdiff(1:ParamNet.CellN,ParamNet.GroupTargetCell);
% IncludeInd=intersect(IncludeInd,find(PosRate>0.01));
% ResParam.yLim=[0 ResponseLim(2)];
% 
% [Res,r,p]=LuPairRegressPlot(ScoreCell(IncludeInd),PosRate(IncludeInd)',ResParam);
% 
% hold on;
% 
% IncludeInd=setdiff(1:ParamNet.CellN,ParamNet.GroupTargetCell);
% IncludeInd=intersect(IncludeInd,find(NegRate<-0.01));
% ResParam.yLim=[ResponseLim(1) 0];
% [Res(2),r(2),p(2)]=LuPairRegressPlot(ScoreCell(IncludeInd),NegRate(IncludeInd)',ResParam);
% 

set(gca,'ylim',ResponseLim,'xlim',ScoreLim);


tHandels={t1 t2 t3 t4};
   % ResParam.Color=[0.1 0.1 0.1];
   % ResParam.yLabel='Inhibition';
   % ResParam.yLim=[ResponseLim(1) 0];
   % 
% t5=subplotPosLu(xMat,yMat,1,5,P);
% [Res(2),r(2),p(2)]=LuPairRegressPlot(ScoreCell(:),NegRate(:),ResParam);
% 
% set(gca,'ytick',ResParam.yLim,'xtick',union(0,ResParam.xLim))

end


function theta=ClockWiseGraph1(G,p,radius)

% Get the number of nodes
N = size(G.Nodes, 1)-1;

% Define the circle radius
% radius = 10; % Adjust as needed

% Compute angles (starting from 12 o'clock position)
angles = linspace(pi/2, -3*pi/2, N+1);  % Starts at 12 o'clock and goes clockwise

% Compute new X, Y coordinates
X_new = radius * cos(angles(1:N));
Y_new = radius * sin(angles(1:N));

theta=angles(1:N);

% Update graph layout
p.XData = [X_new 0];
p.YData = [Y_new 0];

end
