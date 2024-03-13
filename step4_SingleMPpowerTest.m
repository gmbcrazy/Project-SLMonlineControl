
ProcessFolder='F:\LuSLMOnlineTest\03082024\'
RandomDelayInterval=[0 1]; %%Random delay is induced after each trial of stimulation.
PointRepetition=2;  %%Trial Number per each xml MarkPoint stimulation.
nPlane=3;
frameRepetion=40; %%Total repepitions of Z series in T series setting.
MaxFrame=nPlane*frameRepetion;
PV_LinkExcuteFolder(ProcessFolder,RandomDelayInterval,PointRepetition,MaxFrame);