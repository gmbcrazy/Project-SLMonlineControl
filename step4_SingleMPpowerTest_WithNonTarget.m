clear all
ProcessFolder='F:\LuSLMOnlineTest\04222024\SingleP\30PixelFromEdgeExc\';
RandomDelayInterval=[0 1]; %%Random delay is induced after each trial of stimulation.
% PointRepetition=1;  %%Trial Number per each xml MarkPoint stimulation.
nPlane=3;

PreMarkPointRepetition=60;
PostMarkPointRepetition=20;
frameRepetion=PreMarkPointRepetition+PostMarkPointRepetition; %%Total repepitions of Z series in T series setting.
MaxFrame=nPlane*frameRepetion;
BreakPoint=PreMarkPointRepetition*nPlane;
% PV_LinkExcuteFolder(ProcessFolder,RandomDelayInterval,PointRepetition,MaxFrame,BreakPoint);

Round=[3 4 5];
ExcuteIndex=[];
PV_LinkExcuteFolder_byRound(ProcessFolder,RandomDelayInterval,MaxFrame,BreakPoint,Round,ExcuteIndex)


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