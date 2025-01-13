% Load or define your table (TBL)
% Assuming TBL is already loaded in MATLAB workspace
TseriesMultiZFolder='C:\Users\zhangl33\Projects\Project-SLMonlineControl\config\PreGenerateTseriesMultiZ\';
Datalabel='Z11Frame550-10-Jan-2025';

DataFile=[TseriesMultiZFolder Datalabel '.mat'];
load(DataFile);

SaveFolder=[TseriesMultiZFolder Datalabel '\'];
mkdir(SaveFolder)

   P.xLeft=0.1;        %%%%%%Left Margin
   P.xRight=0.02;       %%%%%%Right Margin
   P.yTop=0.1;         %%%%%%Top Margin
   P.yBottom=0.06;      %%%%%%Bottom Margin
   P.xInt=0.02;         %%%%%%Width-interval between subplots
   P.yInt=0.1;         %%%%%%Height-interval between subplots
for iT=1:length(TSeriesBrukerTBL)
   
    figure;
    subplot('position',[0.1 0.35 0.8 0.6]);

TBL = TSeriesBrukerTBL{iT};
% Get unique sequences and functional groups
% Load or define your table (TBL)
% Assuming TBL is already loaded in MATLAB workspace

Xstart=1;
Xwidth=TBL.Reps;
Ystart=0;
Ywidth=1;
Color=[247 150 111;239 109 249;121 247 111]/255;

TimeBarPlot(TBL.SynMPFunGroup,Color,1,Ystart,Xwidth,Ywidth,TBL.VolOut,TBL.PowerZero)

papersizePX=[0 0 20 2];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
saveas(gcf,[SaveFolder 'T' num2str(iT)],'png'); 
% saveas(gcf,[SaveFolder num2str(PSTHparam.PostSLMCal) 'FramePostSLMresponse'],'fig'); 

close all
end




% xlabel('Frames')
