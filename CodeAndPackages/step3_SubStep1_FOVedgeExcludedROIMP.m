% Loading ROI information 
[SavePath,Pos3D,Pos3DRaw,CaData,CaDataPlane,stat,yaml,confSet]=ROIToXYZ(ConfigFolder,ConfigFile);
PlaneZ=confSet.ETL+confSet.scan_Z;
% MultiMatrix3DPlotZ(CaData.PlaneMeanImg,PlaneZ,0.9);
numPlanes=length(confSet.ETL);
iscell=find(CaData.iscell(:,1)==1);


% Intially all cells were dectected by suite2p were considered as SLM targets
SavePathAllPoint=[SavePath 'AllPoint\'];
mkdir(SavePathAllPoint);
IndexNeed=1:1:size(Pos3D,1);
XYZtoMarkPoint(SavePathAllPoint,Pos3D,IndexNeed,yaml,confSet,CaData.statCell);
numPlanes=length(confSet.ETL);
iscell=find(CaData.iscell(:,1)==1);



% Exlude cells near the edge of the FOV as SLM targets
numPoint=size(Pos3D,1);
XYrange=[SLMrangePix;yaml.SLM_Pixels_Y-SLMrangePix];  %%Cell locates close to edge of the view, were not considered as SML targets.
OutRange=find(Pos3D(:,1)<XYrange(1)|Pos3D(:,2)<XYrange(1)|Pos3D(:,1)>XYrange(2)|Pos3D(:,2)>XYrange(2));
CenterCloseI=setdiff(1:numPoint,OutRange);
SavePathExc=[SavePath 'EdgeExc\'];
mkdir(SavePathExc)
% IndexNeed=[1 3 6 13 21 23 28 30 31 32];
% CenterCloseI=CenterCloseI(IndexNeed)
% IndexNeed=1:1:size(Pos3D,1);
XYZtoMarkPoint(SavePathExc,Pos3D,CenterCloseI,yaml,confSet,CaData.statCell);


save([SavePath 'Beh.mat'],'fStim','fSpeed','timeStampCa_Plane','SavePath');