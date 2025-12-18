function [GgraphOut,tHandels]=tblMeta2CircosPlot(tbl_whisktemp, group, GroupParamNet, scoreName)
% plotGroupData - Plot group network and response/score relationships
%
% Syntax:
%   plotGroupData(tbl_whisktemp, group, GroupParamNet, scoreName)
%
% Inputs:
%   tbl_whisktemp : table containing the data
%   group         : numeric or categorical group identifier
%   GroupParamNet : structure containing group parameters (color, limits, maps, etc.)
%   scoreName     : string, name of the score variable in tbl_whisktemp
%
% Example:
%   plotGroupData(tbl_whisktemp, 1, GroupParamNet, 'SpeedR');
ResponseLim = GroupParamNet.ResponseLim;
ScoreLim = GroupParamNet.ScoreLim;
ResponseMap = GroupParamNet.ResponseMap;
ScoreMap = GroupParamNet.ScoreMap;
SessionMap = GroupParamNet.SessionMap;


% Filter table for target group
if group==0    %%If group == 0, use ShamOpto values as reponses.
   tempCelltbl=tbl_whisktemp;
   tempCelltbl(tempCelltbl.TargetCell == 0 & tempCelltbl.NonTargetCell == 0, :) = [];
   tempCelltbl(:,"Response")=[];
   idx = find(tempCelltbl.Properties.VariableNames == "ShamOpto");
   tempCelltbl.Properties.VariableNames{idx} = 'Response';

   tempCelltbl2=groupsummary(tempCelltbl,'Cell', 'mean');
   tempCelltbl2= groupsummaryBack2OldNames(tempCelltbl,tempCelltbl2,'mean');
   tempCelltbl=tempCelltbl2;

else
   tempCelltbl=tbl_whisktemp(tbl_whisktemp.Group==group,:);
   tempCelltbl(tempCelltbl.TargetCell == 0 & tempCelltbl.NonTargetCell == 0, :) = [];
end





% Index cells
CellAll = unique(tempCelltbl.Cell);
cellI = 1:length(CellAll);
tempCelltblIndex = table(cellI(:), CellAll(:), 'VariableNames', {'CellI','Cell'});
tempCelltbl = innerjoin(tempCelltbl, tempCelltblIndex, "Keys", "Cell");

% Sort by score
[~, sortI] = sort(tempCelltbl.(scoreName), 'descend');
tempCelltbl = tempCelltbl(sortI, :);

Adj = NaN(max(cellI)+1);
Adj(end,1:end-1) = tempCelltbl.Response;

Node = [tempCelltbl.Response; 0];
Node(tempCelltbl.TargetCell == group) = 1;

radius = 5;
NodeColor = repmat([0.8 0.8 0.8], max(cellI)+1, 1);
if group~=0
NodeColor(tempCelltbl.TargetCell == group, :) = repmat(GroupParamNet.GroupColor(group,:), sum(tempCelltbl.TargetCell == group), 1);
NodeColor(end,:) = GroupParamNet.GroupColor(group,:);
end


Adj(end, tempCelltbl.TargetCell ~= group) = NaN;
Adj(end, tempCelltbl.NonTargetCell == 1) = NaN;

figure;
CirSubP = GroupParamNet.CirSubP;

% Network plot
t1 = subplotLU(1,2,1,1,CirSubP);
[GgraphOut.G, GgraphOut.p] = MarkovState_HeatStrPlot(abs(Adj), zeros(size(Adj,1),1)+0.1, Adj, Node, NodeColor, ResponseMap, ResponseLim);
theta = ClockWiseGraph1(GgraphOut.G, GgraphOut.p, radius);

% Circos bands
clear CircosBand;
CircosBand(1).rband = [radius+0.5; radius+2];
CircosBand(1).Values = tempCelltbl.Response;
CircosBand(1).Clim = ResponseLim;
CircosBand(1).theta = theta;
CircosBand(1).Colormap = ResponseMap;
CircosBand(1).arc_gap = (theta(2)-theta(1))/10;

CircosBand(end+1).rband = CircosBand(end).rband + 1.6;
CircosBand(end).Values = tempCelltbl.(scoreName);
CircosBand(end).Clim = ScoreLim;
CircosBand(end).theta = theta;
CircosBand(end).Colormap = ScoreMap;
CircosBand(end).arc_gap = 0;

CircosBand(end+1).rband = CircosBand(end).rband + 1.6;
CircosBand(end).rband(2) = CircosBand(end).rband(1) + 0.4;
CircosBand(end).Values = tempCelltbl.Session;
CircosBand(end).Clim = [min(tempCelltbl.Session);max(tempCelltbl.Session)];
CircosBand(end).theta = theta;
CircosBand(end).Colormap = SessionMap;
CircosBand(end).arc_gap = 0;

AddCircosBand(CircosBand);
GgraphOut.p.ArrowSize = 10;
GgraphOut.p.MarkerSize = repmat(0.5, size(NodeColor,1), 1);
GgraphOut.p.MarkerSize(tempCelltbl.TargetCell == group) = 3;
GgraphOut.p.MarkerSize(end) = 10;
GgraphOut.p.NodeColor = NodeColor;

colormap(t1, ScoreMap);
set(gca, 'Clim', ScoreLim);
bar1 = colorbar(t1);
bar1.Location = 'northoutside';
bar1.Position = [t1.Position(1)+0.1, t1.Position(2)+t1.Position(4)+0.01, 0.1, 0.02];
bar1.Label.String = scoreName; % Updated to use input variable name
bar1.Ticks = union(ScoreLim, 0);

% Regression plot
t2 = subplotLU(1,2,1,2,CirSubP);



ParamTemp.Marker='o';
ParamTemp.MarkerSize=8;
ParamTemp.Rtype='pearson';
% ParamTemp.xLim=[min(data1) max(data1)];
% ParamTemp.yLim=[min(data2) max(data2)];
ParamTemp.xLabel=scoreName;
ParamTemp.yLabel='Response residuals';
ParamTemp.RawcorrPlot = 0;
ParamTemp.xLim=[];
ParamTemp.yLim=[];


ParamTemp.Color = Value2Color(unique(tempCelltbl.Session), SessionMap, [min(tempCelltbl.Session);max(tempCelltbl.Session)]);
tblCovPlot=tempCelltbl;
tblCovPlot=tblCovPlot(tblCovPlot.NonTargetCell==1,:);
LuPairRegressPlot_Group_Cov(tblCovPlot.(scoreName), tblCovPlot.Response, tblCovPlot.Speed, tblCovPlot.Session, ParamTemp);
colormap(t2, ResponseMap);
set(gca, 'clim', ResponseLim);
bar2 = colorbar(t2);
bar2.Location = 'northoutside';
bar2.Position = [t1.Position(1)+0.25, t1.Position(2)+t1.Position(4)+0.01, 0.1, 0.02];
bar2.Label.String = 'Response';
bar2.Ticks = union(ResponseLim, 0);

% Session colorbar
t3 = subplot('Position', [t2.Position(1)+t2.Position(3)+0.02, 0.4, 0.02, 0.3]);
set(t3, 'Visible', 'off')

colormap(t3, ParamTemp.Color);
bar3 = colorbar(t3);
bar3.Location = 'eastoutside';
bar3.Position = [t3.Position(1)+0.0, t3.Position(2), 0.02, t3.Position(4)];

bar3.Label.String = 'Session';
bar3.Ticks = [];

tHandels={t1 t2 t3};
end
