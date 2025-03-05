
clear all
Param.TrialRepTh=4;
Param.PreImgN=10;
% PostImgN=min(90,median(DataList.totalRepetitions));

Param.PreSLMFrameN=20; %Frame num. before SLM to calculate baseline level
Param.PostSLMFrameN=3; %Frame num. after SLM to calculate response level
Param.TrialRepTh=4; %%Minimal num. of trials to quantify SLM response
Param.DistTh=10;%%Maximal distance by pixels considered between target cell and SLM target center.

TargetZ=[];
SessinInfo=[];
GlobalSavePath='\\nimhlabstore1.nimh.nih.gov\UFNC\FNC2\Zhang\SLM\SummarySLMResponse\AllCurentSLMtarget';


ProcessFolder='F:\LuSLMOnlineTest\SL0543-Ai203\09102024\SingleP\Top8SpeedStimEdgeExc\';
SessinInfo=[];
SLMTable=SLMtargetAnalysis(ProcessFolder, Param,GlobalSavePath,SessinInfo, TargetZ);


ProcessFolder='F:\LuSLMOnlineTest\SL0242-Ai203\09042024\SingleP\Top5SpeedStimEdgeExc\';
SessinInfo=[];
SLMTable=SLMtargetAnalysis(ProcessFolder, Param,GlobalSavePath,SessinInfo, TargetZ);


ProcessFolder='F:\LuSLMOnlineTest\SL0242-Ai203\08292024\SingleP\Top5SpeedStimEdgeExc\';
SessinInfo=[];
SLMTable=SLMtargetAnalysis(ProcessFolder, Param,GlobalSavePath,SessinInfo, TargetZ);



ProcessFolder='D:\LuSLMOnlineTest\SL0340\05302024\Ch2\TSeries\PowerLevel1\';
SessinInfo=1;
SLMTable=SLMtargetAnalysis(ProcessFolder, Param,GlobalSavePath,SessinInfo, TargetZ);


ProcessFolder='\\nimhlabstore1.nimh.nih.gov\UFNC\FNC3Z\Sutie2p-Processed\SL1865-Ai203\05142024\SingleP\20PixelFromEdgeExc\';
SessinInfo=1;
SLMTable=SLMtargetAnalysis(ProcessFolder, Param,GlobalSavePath,SessinInfo, TargetZ);



ProcessFolder='F:\LuSLMOnlineTest\SL0541-Emx1Ai96\09122024\SingleP\Top8SpeedStimEdgeExc\';
SessinInfo=[];
SLMTable=SLMtargetAnalysis(ProcessFolder, Param,GlobalSavePath,SessinInfo, TargetZ);

ProcessFolder='F:\LuSLMOnlineTest\SL0541-Emx1Ai96\09172024\SingleP\Top10SpeedStimEdgeExcTest\';
SessinInfo=[];
SLMTable=SLMtargetAnalysis(ProcessFolder, Param,GlobalSavePath,SessinInfo, TargetZ);



ProcessFolder='F:\LuSLMOnlineTest\SL0702-G7-Ai203\09192024\SingleP\Top15SpeedStimEdgeExc\';
SessinInfo=[];
SLMTable=SLMtargetAnalysis(ProcessFolder, Param,GlobalSavePath,SessinInfo, TargetZ);
