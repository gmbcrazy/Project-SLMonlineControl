% clear all
SaveFunCon='\\nimhlabstore1.nimh.nih.gov\UFNC\FNC2\Zhang\Projects\Project-LocalProcessing\Step4\awakeRefSpon\GroupSLM10Sessions\';
T=readtable([SaveFunCon 'SLMGroupResTrialDynWin3.csv']);




% SLMGroupTableTrial = readtable('SLMGroupResTrialDynWin1.csv');

% Tavg = varfun(@mean, SLMGroupTableTrial, ...
%                'InputVariables', 'Response', ...
%                'GroupingVariables', {'Session','Cell'});

% Tavg.Properties.VariableNames{'mean_Response'} = 'AvgResponse';
% Tavg now has Session, Cell, and AvgResponse

% 1) Read your table

% 2) Extract the one‐row “key” table of constant vars
keyVars = {'Session','Cell', ...
           'SpeedR','StimR','TargetCell', ...
           'Group','Whisk','PowerZero','AwakeState','TargetSpeedR','TargetStimR','Speed'};  
KeyTab = unique( T(:, keyVars) , 'rows' );

KeyPar={'Whisk','PowerZero','AwakeState','TargetCell'};

KeyPar = unique( T(:, keyVars) , 'rows' );


% 3) Compute average Response by (Session,Cell)
MeanTab = groupsummary( T, {'Session','Cell','Group','Whisk','PowerZero'}, 'mean', 'Response');
MeanTab.Properties.VariableNames{'mean_Response'} = 'AvgResponse';
MeanTab.GroupCount = [];   % drop the count column

% 4) Join them back so every row has AvgResponse plus your other columns
AveTable = join( KeyTab, MeanTab, 'Keys', {'Session','Cell','Whisk','Group','PowerZero'});

Session=unique(AveTable.Session);



NonWhiskT=AveTable(AveTable.Whisk==1&AveTable.PowerZero==0,:);
NonTargetY=groupsummary(NonWhiskT(NonWhiskT.TargetCell==0,:),{'Session','Group'}, @(x,y) corr(x,y),{"AvgResponse","SpeedR"});

% A=NonWhiskT(NonWhiskT.TargetCell==1,:);
TargetX=groupsummary(NonWhiskT(NonWhiskT.TargetCell==1,:),{'Session','Group'}, 'mean',"TargetSpeedR");
Cov=groupsummary(NonWhiskT(NonWhiskT.TargetCell==1,:),{'Session','Group'}, 'mean',"Speed");


NonTargetY=groupsummary(NonWhiskT(NonWhiskT.TargetCell==0,:),{'Session','Group'}, @(x,y) corr(x,y),{"AvgResponse","StimR"});

% A=NonWhiskT(NonWhiskT.TargetCell==1,:);
TargetX=groupsummary(NonWhiskT(NonWhiskT.TargetCell==1,:),{'Session','Group'}, 'mean',"TargetStimR");
Cov=groupsummary(NonWhiskT(NonWhiskT.TargetCell==1,:),{'Session','Group'}, 'mean',"Speed");

% B=groupsummary(A,{'Session','Group'}, 'mean','TargetSpeedR');


% A=NonWhiskT(NonWhiskT.TargetCell==1,:);

ProcessPar.PlotFigure=0;
ProcessPar.GroupLabel={'L','S','N'};
ProcessPar.GroupList=[1 2 3];
ProcessPar.GroupColor=[255 51 153;91 20 212;121 247 111]/255;
ProcessPar.PowerZeroColor=[0.5 0.5 0.5];
ProcessPar.PowerZero=[0 1];
ProcessPar.PowerZeroLabel={'SLM','FakeSLM'};
ProcessPar.VolOut=[0 1];
ProcessPar.VolOutLabel={'NoSensory','Sensory'}
ProcessPar.AwakeState=[1];
ProcessPar.AwakeStateLabel={'Awake'};
% % GroupMetaName=[ProcessPar.GroupLabel {'FakeSLM'}];

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
   Param.xLim=[-0.2 0.5];
   Param.yLim=[-0.2 0.5];
   

NonWhiskT(NonWhiskT.TargetCell==1,:)


% [OutPut,r,p]=LuPairRegressPlot_Group(TargetX.mean_TargetSpeedR,NonTargetY.fun1_AvgResponse_SpeedR,TargetX.Group,Param,Cov.mean_Speed)    
% figure;
% [OutPut,r,p]=LuPairRegressPlot_Group_Cov(TargetX.mean_TargetSpeedR,NonTargetY.fun1_AvgResponse_SpeedR,Cov.mean_Speed,TargetX.Group,Param)    



   Param.xLim=[-0.3 0.6];
   Param.yLim=[-0.3 0.6];
  
figure;
[OutPut,r,p,h]=LuPairRegressPlot_Group_Cov(TargetX.mean_TargetStimR,NonTargetY.fun1_AvgResponse_StimR,Cov.mean_Speed,TargetX.Group,Param)    


% Data1=groupsummary(NonWhiskT,{'Session','Group'}, @(x,y) corr(x,y),{"mean_Response","SpeedR"});





% ----
% Result now has columns:
%   Session | Cell | SpeedR | StimR | TargetCell | Group | Whisk | PowerZero | AwakeState | AvgResponse
