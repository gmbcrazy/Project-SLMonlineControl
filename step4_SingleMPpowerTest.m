clear all
ProcessFolder='F:\LuSLMOnlineTest\04092024\'
RandomDelayInterval=[0 1]; %%Random delay is induced after each trial of stimulation.
PointRepetition=2;  %%Trial Number per each xml MarkPoint stimulation.
nPlane=3;

PreMarkPointRepetition=25;
PostMarkPointRepetition=30;
frameRepetion=PreMarkPointRepetition+PostMarkPointRepetition; %%Total repepitions of Z series in T series setting.
MaxFrame=nPlane*frameRepetion;
BreakPoint=PreMarkPointRepetition*nPlane;
PV_LinkExcuteFolder(ProcessFolder,RandomDelayInterval,PointRepetition,MaxFrame,BreakPoint);


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