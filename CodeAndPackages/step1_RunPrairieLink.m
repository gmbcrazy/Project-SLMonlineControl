ConfigFolder='C:\Users\User\Project-SLMonlineControl\config\';


%% Test whisk stim
disp('turn off PMT and power before test')
TSeriesENVFile=dir([ConfigFolder 'WhiskerStimImgSeqTest.env'])
LoadTSeriestoBruker(TSeriesENVFile)


%% Get motion correction refImg
TSeriesENVFile=dir([ConfigFolder 'RefZforMotionCorrection.env'])
LoadTSeriestoBruker(TSeriesENVFile)



%% Spontnous beh recording
TSeriesENVFile=dir([ConfigFolder 'SponBehvZ3F4100.env'])
LoadTSeriestoBruker(TSeriesENVFile)
PrairieLink_RawDataStream_UFNC