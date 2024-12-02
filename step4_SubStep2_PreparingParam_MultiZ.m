
nPlane=length(confSet.ETL);
frameRepetion=PreMarkPointRepetition+PostMarkPointRepetition; %%Total repepitions of Z series in T series setting.


PVparam.maxFrame=nPlane*frameRepetion;
PVparam.BreakPointFrame=PreMarkPointRepetition*nPlane;
PVparam.maxFrame=nPlane*frameRepetion;

%param for calculate the PSTH heatmap for online analysis
PSTHparam.PreInd=PreMarkPointRepetition-PreSLMCal:PreMarkPointRepetition;
PSTHparam.PostInd=PreMarkPointRepetition+1:PreMarkPointRepetition+PostSLMCal;
PSTHparam.Plot=1;
PSTHparam.SmoothSD=1;
PSTHparam.ColorMap=colorMapPN1;
PSTHparam.Clim=[-400 400];

%param for xml files
XMLparam.ProcessFolder=ProcessFolder;
XMLparam.TotalRounds=confSet.NumTrial;
PointAll=1:size(Pos3Dneed,1);

XMLparam.Laser=1.5;               %<-------------------------------------------------------------------------------Edit, starting laser power to test    
XMLparam.RoundID=1;               %starting round
PointAll=1:size(Pos3Dneed,1);     %All possible MarkPoints for testing
PointsTest=PointAll;              %Initial test Points, this would be updated automatically later


%param for ROI neighbourhood to determine wether there is SLM response.
ROIparam.TotalSLMPos3D=Pos3Dneed;    %%such that ROIparam.PointAll=1:size(Pos3Dneed,1)
ROIparam.PointAll=PointAll;
ROIparam.PlaneZ=PlaneZ;
ROIparam.CellSize=20;                %%normal neuron diameter by um;        
ROIparam.threshold_percentage=0.3;   %%thereshold to define responsive fields SLM responsive heatmap: percentage*Peak rate
ROIparam.thNum=15;                   %%Minimal single responsive field by pixels
ROIparam.max_distance=ceil(ROIparam.CellSize*2/3/umPerPixel);  %% 2/3 diameter of a cell by pixel as maximal response region-SLM center distance
ROIparam.min_region_size=10;
ROIparam.PeakTh=200;
ROIparam.min_merged_region_size=20;  %%Minimal total size of responsive fields by pixels
ROIparam.contourMethod='perim';      %%Method to detect ROI boader of responsive fields
ROIparam.NeighbourHfWidthPixel=20;   %%PixFromMedCenter: Number of pixels from the median center of each ROI to get the ROI neighborhood.
ROIparam.umPerPixel=mean([yaml.umPerlPixelX yaml.umPerlPixelY]);
ROIparam.Colormap=colorMapPN1;                  
% ROIparam.LaserPower=confSet.UncagingLaserPower;   


