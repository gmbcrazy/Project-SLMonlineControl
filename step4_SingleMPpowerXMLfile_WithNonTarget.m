clear all
% ProcessFolder='F:\LuSLMOnlineTest\04222024\SingleP\30PixelFromEdgeExc\';
% load('C:\Users\zhangl33\Projects\GenMatCode\Plotfun\Color\colorMapPN3.mat');

load('C:\Users\User\Project-SLMonlineControl\subfun\Color\colorMapPN3.mat');

% ProcessFolder='F:\LuSLMOnlineTest\04222024\SingleP\30PixelFromEdgeExc\';
ProcessFolder='F:\LuSLMOnlineTest\MouseMei03\04242024Test\50PixelFromEdgeExc\';

DataFolder=[ProcessFolder 'Data\'];

% load([ProcessFolder 'SLMIncludedIndFromIscell.mat'])

% RandomDelayInterval=[0 1]; %%Random delay is induced after each trial of stimulation.
% PointRepetition=1;  %%Trial Number per each xml MarkPoint stimulation.
nPlane=length(confSet.ETL);

PreMarkPointRepetition=60;
PostMarkPointRepetition=20;
frameRepetion=PreMarkPointRepetition+PostMarkPointRepetition; %%Total repepitions of Z series in T series setting.
PVparam.maxFrame=nPlane*frameRepetion;
PVparam.BreakPointFrame=PreMarkPointRepetition*nPlane;


XMLparam.Point=2;
XMLparam.Laser=1.5;
XMLparam.RoundID=1;
XMLparam.ProcessFolder=ProcessFolder;
load([XMLparam.ProcessFolder 'SLMIncludedIndFromIscell.mat']);
PSTHparam.TargetPos=Pos3DNeed(XMLparam.Point,:);
PSTHparam.PreInd=PreMarkPointRepetition-20:PreMarkPointRepetition;
PSTHparam.PostInd=PreMarkPointRepetition+1:PreMarkPointRepetition+3;
PSTHparam.Plot=1;
PSTHparam.SmoothSD=1;
PSTHparam.ColorMap=ColorPN3;
PSTHparam.Clim=[-200 200]

% PV_LinkExcuteFolder(ProcessFolder,RandomDelayInterval,PointRepetition,MaxFrame,BreakPoint);

Round=[1];
% ExcuteIndex=[];



% PV_LinkExcuteXML(ProcessFolder,RandomDelayInterval)
PV_LinkExcuteXML(XMLparam,PVparam,confSet)
PV_LinkExcuteXML(XMLparam,PVparam,confSet,PSTHparam);

% xmlCSV=[ProcessFolder 'xmlListCurrent.csv'];
% PV_LinkExcuteFolder_byRound(xmlCSV,RandomDelayInterval,MaxFrame,BreakPoint,Round,[2]);

%Point 7 Only

%%disable breakP, 2 files?


%%Enable breakP, 2files.  Flush loop.


%Add Saline
%%Flush once after breakP.  Flush once.
%File 11





%Add Saline
%%Flush once at the break point;
%File 13

%File 17



%Point 6 Only
%%Flush once at the break point;
%File 13

%File 17