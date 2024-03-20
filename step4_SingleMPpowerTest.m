
ProcessFolder='F:\LuSLMOnlineTest\03202024\'
RandomDelayInterval=[0 1]; %%Random delay is induced after each trial of stimulation.
PointRepetition=2;  %%Trial Number per each xml MarkPoint stimulation.
nPlane=3;
MaxFrame=nPlane*frameRepetion;

PreMarkPointRepetition=11;
PostMarkPointRepetition=16;
frameRepetion=PreMarkPointRepetition+PostMarkPointRepetition; %%Total repepitions of Z series in T series setting.
MaxFrame=nPlane*frameRepetion;
BreakPoint=PreMarkPointRepetition*nPlane;
PV_LinkExcuteFolder(ProcessFolder,RandomDelayInterval,PointRepetition,MaxFrame,BreakPoint);