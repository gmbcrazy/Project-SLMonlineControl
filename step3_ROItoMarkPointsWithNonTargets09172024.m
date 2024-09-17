clear all
ConfigFolder='C:\Users\User\Project-SLMonlineControl\config\';
% ConfigFolder='C:\Users\zhangl33\Projects\Project-SLMonlineControl\config\';
[SavePath,Pos3D,Pos3DRaw,CaData,CaDataPlane,stat,yaml,confSet]=ROIToXYZ(ConfigFolder);

PlaneZ=confSet.ETL+confSet.scan_Z;
% MultiMatrix3DPlotZ(CaData.PlaneMeanImg,PlaneZ,0.9);
numPlanes=length(confSet.ETL);
iscell=find(CaData.iscell(:,1)==1);

if ~exist('fSpeed')
    [fSpeed,fStim,timeStampCa_Plane]=PV_VolExtract(confSet);
    fSpeed=AmpNormalizeRow(double(fSpeed)',[0 100])';
    fStim=AmpNormalizeRow(double(fStim)',[0 100])';
end


XTimesStdTh = 3;
MinInterVal = 20;
maxLag = 10;


% Intially all cells were dectected by suite2p were considered as SLM targets
SavePathAllPoint=[SavePath 'AllPoint\']
mkdir(SavePathAllPoint)
IndexNeed=1:1:size(Pos3D,1);
XYZtoMarkPoint(SavePathAllPoint,Pos3D,IndexNeed,yaml,confSet,CaData.statCell);
numPlanes=length(confSet.ETL);
iscell=find(CaData.iscell(:,1)==1);



%% Exlude cells near the edge of the FOV as SLM targets
SLMrangePix=20; %Pixel number close to FOV is excluded
numPoint=size(Pos3D,1);
XYrange=[SLMrangePix;yaml.SLM_Pixels_Y-SLMrangePix]  %%Cell locates close to edge of the view, were not considered as SML targets.
OutRange=find(Pos3D(:,1)<XYrange(1)|Pos3D(:,2)<XYrange(1)|Pos3D(:,1)>XYrange(2)|Pos3D(:,2)>XYrange(2));
CenterCloseI=setdiff(1:numPoint,OutRange);
SavePathExc=[SavePath 'EdgeExc\']
mkdir(SavePathExc)
% IndexNeed=[1 3 6 13 21 23 28 30 31 32];
% CenterCloseI=CenterCloseI(IndexNeed)
% IndexNeed=1:1:size(Pos3D,1);
XYZtoMarkPoint(SavePathExc,Pos3D,CenterCloseI,yaml,confSet,CaData.statCell);


%% Automatic generate Non-Targets 
NonTargetPath=[SavePath  'NonTargetsTest\'];
mkdir(NonTargetPath)
SaveNonTargets=[NonTargetPath 'Raw']
BloodVesselTh=0.1;
[NonTargets,NonTargetsPlane]=NonTargetGeneration(SaveNonTargets,Pos3DRaw,CaData.CellPlaneIDRaw,yaml,confSet,CaData.PlaneMeanImg,BloodVesselTh);
% confSet.NumNonTarget=20;
% confSet.RadiusAvoidParam=3;
% % [NonTargets,NonTargetsPlane]=NonTargetGeneration(SaveNonTargets,Pos3DRaw,CaData.CellPlaneIDRaw,yaml,confSet);
figure;
FigSavePath=[NonTargetPath 'Raw'];
PlotTargetNonTarget(Pos3D,NonTargets,NonTargetsPlane,CaData,FigSavePath)


% [ s ] = gpl2struct('F:\LuSLMOnlineTest\04112024\SingleP\NonTargets\NonTargets.gpl')

%% Please do mannual correction to exclude disqualified non-targets in PV, after that, exported all selected targets to SelectedFromRaw.gpl file
s=gpl2struct([NonTargetPath 'SelectedFromRaw.gpl']);
% TempTable=struct2table(s.PVGalvoPointList.PVGalvoPoint{1}.Attributes);
for i=1:length(s.PVGalvoPointList.PVGalvoPoint)
    TempStr(i)=s.PVGalvoPointList.PVGalvoPoint{i}.Attributes;
end
NonTargetNeed = convertTableEntries(struct2table(TempStr));
NonTargetNeedInd = NonTargetNeed.Index+1;

NonTargets=NonTargets(NonTargetNeedInd,:);
NonTargetsPlane=NonTargetsPlane(NonTargetNeedInd,:);
SaveNonTargets=[NonTargetPath 'SelectedNonTargets'];
% [NonTargets,NonTargetsPlane]=NonTargetGeneration(SaveNonTargets,NonTargets,NonTargetsPlane,yaml,confSet);
MarkPoints3D_GPLmaker(NonTargets, yaml, true, confSet.SpiralSizeUM, confSet.SpiralRevolution,SaveNonTargets,[], 'NonTarget');    %%<< This file is not used in this script, just saved for note, or visulization in PV.

% load([SavePath 'NonTargets\SelectedNonTargets']);

% Visulize all targets and non-targets
figure;
FigSavePath=[NonTargetPath 'TargetNonTarget'];
PlotTargetNonTarget(Pos3D,NonTargets,NonTargetsPlane,CaData,FigSavePath)



% Visulize all targets and non-targets selected for each trial 
for iTrial=1:confSet.NumTrial
    figure;
    IndexNonTargetTrial(:,iTrial)=randperm(size(NonTargets,1),confSet.NumNTperTrial);
    FigSavePath=[NonTargetPath 'Trial' num2str(iTrial)];
    PlotTargetNonTarget(Pos3D,NonTargets(IndexNonTargetTrial(:,iTrial),:),NonTargetsPlane(IndexNonTargetTrial(:,iTrial)),CaData,FigSavePath);
    close all
end





% % %% Intially all cells were dectected by suite2p were considered as SLM targets
% SavePathAllPoint=[SavePath 'AllPoint\']
% mkdir(SavePathAllPoint)
% IndexNeed=1:1:size(Pos3D,1);
% % XYZtoMarkPoint(SavePathAllPoint,Pos3D,IndexNeed,yaml,confSet);
% XYZtoMarkPoint_NT(SavePathAllPoint,Pos3D,IndexNeed,NonTargets, IndexNonTargetTrial, yaml,confSet)
% 




%% Exlude cells near the edge of the FOV as SLM targets
% SLMrangePix=50; %Pixel number close to FOV is excluded
% numPoint=size(Pos3D,1);
% XYrange=[SLMrangePix;yaml.SLM_Pixels_Y-SLMrangePix]  %%Cell locates close to edge of the view, were not considered as SML targets.
% OutRange=find(Pos3D(:,1)<XYrange(1)|Pos3D(:,2)<XYrange(1)|Pos3D(:,1)>XYrange(2)|Pos3D(:,2)>XYrange(2));
% CenterCloseI=setdiff(1:numPoint,OutRange);
SavePathExc=[SavePath num2str(SLMrangePix) 'PixelFromEdgeExc\']
mkdir(SavePathExc)
XYZtoMarkPoint_NT_PairGplXml(SavePathExc, Pos3D, CenterCloseI, NonTargets, IndexNonTargetTrial, yaml, confSet,CaData.statCell);

%% 

% [Neighbourhood1,  ~] = MarkPoint2Neighbourhood(MarkPoints, radius*2, numPlanes,planeSize);
% [Neighbourhood2,  ~] = MarkPoint2Neighbourhood(newMarkPoints, radius, numPlanes,planeSize);
% 




%% Including top cells highly correlated associated with speed.

deltaFoF=F2deltaFoF(double(CaData.F),double(CaData.Fneu),double(confSet.fs));
NeuroData=AmpNormalizeRow(deltaFoF',[0 100])';

for iPlane=1:numPlanes
    I1=find(CaData.CellPlaneID==iPlane);
    [rSpeed(I1),pSpeed(I1)]=corr(NeuroData(:,iscell(I1)),fSpeed(:,iPlane),'type','Spearman','rows','pairwise');
end
clear rStim
for iCell = 1:size(iscell, 1)
    iPlane=CaData.CellPlaneID(iCell);
    [c, lags] = xcorr(NeuroData(:,iscell(iCell)), fStim(:,iPlane), maxLag, 'coeff');
    PostI = find(lags >= 0);
    [~, i1] = max(abs(c(PostI)));
    rStim(iCell, 1) = c(PostI(i1));
end
[~,rCenterIStim]=sort(rStim(CenterCloseI),'descend');  %%Cell locates close to edge of the view, were not considered as SML targets.
[~,~,rankStim]=intersect(1:length(CenterCloseI),rCenterIStim);

% ExcludeStimI=setdiff(CenterCloseI,)
[~,rCenterISpeed]=sort(rSpeed(CenterCloseI),'descend');  %%Cell locates close to edge of the view, were not considered as SML targets.
[~,~,rankSpeed]=intersect(1:length(CenterCloseI),rCenterISpeed);


TopCellN=6;  
IncludeCellI=union(rCenterISpeed(1:TopCellN),rCenterIStim(1:TopCellN));
rScore=double([rSpeed(:) rStim(:)]);
rScore=rScore(CenterCloseI,:);
rRank=[rankSpeed(:) rankStim(:)]
rScoreInclude=rScore(IncludeCellI,:);
rRankInclude=rRank(IncludeCellI,:);

IncludeCellFunFilter=CenterCloseI(IncludeCellI);
SavePathStimSpeed=[SavePath 'Top' num2str(TopCellN) 'SpeedStimEdgeExcTest\']
mkdir(SavePathStimSpeed)
% XYZtoMarkPoint_NT(SavePathSpeed,Pos3D,IncludeCellFunFilter,NonTargets, IndexNonTargetTrial, yaml,confSet)
% XYZtoMarkPoint_NT(SavePathSpeed,Pos3D,IncludeCellFunFilter,NonTargets, IndexNonTargetTrial, yaml,confSet)
IncludePath=[SavePathStimSpeed '\AllIncludedOrigin\'];
mkdir(IncludePath)
XYZtoMarkPoint(IncludePath,Pos3D,IncludeCellFunFilter,yaml,confSet,CaData.statCell);


EditedTable=gpl2Table([IncludePath 'SelectedGPL.gpl'])
OriTable=gpl2Table([IncludePath 'GPL.gpl'])
[temp,I1,I2]=intersect(EditedTable.Name,OriTable.Name);
IncludePath=[SavePathStimSpeed '\AllIncluded\'];
mkdir(IncludePath)
IncludeCellFinal=IncludeCellFunFilter(I2);
SLMIncludedIndFromIscell=IncludeCellFinal;
[XYPosPixel,Z]=gplXYtoPixel(EditedTable,yaml);
Pos3DNeed=[XYPosPixel Z(:,1)];
Cellstat=CaData.statCell;
save([IncludePath 'SLMIncludedIndFromIscell.mat'],'SLMIncludedIndFromIscell','Pos3DNeed','yaml','confSet','NonTargets','IndexNonTargetTrial','Cellstat');

Pos3DFromGPL=[EditedTable.X EditedTable.Y EditedTable.Z];

NonTargetsTable=gpl2Table([NonTargetPath 'SelectedFromRaw.gpl']);
NonTargetsGPL=[NonTargetsTable.X NonTargetsTable.Y NonTargetsTable.Z];


SavePathStimSpeed=[SavePath 'FinalTop' num2str(TopCellN) 'SpeedStimEdgeExc\'];
mkdir(SavePathStimSpeed);
% XYZtoMarkPoint_NT_PairGplXml(SavePathStimSpeed, Pos3D, IncludeCellFunFilter, NonTargets, IndexNonTargetTrial, yaml, confSet,CaData.statCell);
EditedGPLtoMarkPoint_NT_PairGplXml(SavePathStimSpeed, Pos3DFromGPL, SLMIncludedIndFromIscell, NonTargetsGPL, IndexNonTargetTrial, yaml, confSet,Cellstat)



%%
figure;
plot(rSpeed(CenterCloseI),'r-');
hold on;
plot(rStim(CenterCloseI),'g-');
plot(rCenterISpeed(1:TopCellN),0.3,'ro')
plot(rCenterIStim(1:TopCellN),0.3,'g*')



figure;
for iTop=1:TopCellN
% subplot(TopCellN,1,i)
tempData=NeuroData(:,iscell(CenterCloseI(rCenterIStim(iTop))));
tempData=SmoothDec(tempData,1);
plot(tempData+iTop);
hold on;
plot(mean(fStim,2)*(TopCellN+1),'r');
end



figure;
for iTop=1:TopCellN
% subplot(TopCellN,1,i)
tempData=NeuroData(:,iscell(CenterCloseI(rCenterISpeed(iTop))));
tempData=SmoothDec(tempData,1);
plot(tempData+iTop);
hold on;
plot(mean(fSpeed,2)*(TopCellN+1),'r');
end


subplot(2,1,1)
imagesc(NeuroData(:,iscell(CenterCloseI))');
subplot(2,1,2)
plot(mean(fSpeed,2));hold on;


figure;
for i=1:length(CenterCloseI)
    subplot(length(CenterCloseI),1,i)
plot(NeuroData(:,iscell(CenterCloseI(i))));hold on;
plot(mean(fStim,2)+5,'r-');
end


figure;
plot(rSpeed,'r');
hold on;
plot(rStim,'g');
hold on;
plot(rFspks,'b')

figure;
plot(-log10(pFspks),'b')
[~,r1]=sort(rFspks(CenterCloseI),'descend')  %%Cell locates close to edge of the view, were not considered as SML targets.
TopCellN=8;         %%   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<Top X cells highly associated with speed 
TopSpeedCellI=CenterCloseI(r1(1:TopCellN));
SavePathSpeed=[SavePath 'Top' num2str(TopCellN) 'SpeedEdgeExc\']
mkdir(SavePathSpeed)
XYZtoMarkPoint_NT(SavePathSpeed,Pos3D,TopSpeedCellI,NonTargets, IndexNonTargetTrial, yaml,confSet)





%% 
% load('C:\Users\zhangl33\Projects\GenMatCode\Plotfun\Color\colorMapPN3.mat');
% rCellspks=partialcorr(CaData.spks(iscell,:)',fSpeed);
rCellspks=corr(double(CaData.spks(iscell,:))','type','Spearman','rows','pairwise');
load('C:\Users\zhangl33\Projects\GenMatCode\Plotfun\Color\colorMapPN3.mat')
AdjMat=rCellspks+1;  %%%%%%%%%%Non-negative weighted correlation.
clear SampleCorr;
for itt=1:size(AdjMat,1)
    AdjMat(itt,itt)=0;%%%%%%diagonal zeros for Adjcent matrix
end
k = full(sum(AdjMat));
twom = sum(k);
B = full(AdjMat - gamma*k'*k/twom);
tic
k = full(sum(AdjMat));
twom = sum(k);
gamma = 1.05;  %%Community Clustering Paramter%%
limit =100000; %%memory consideration for community clustering 
B = @(i) AdjMat(:,i) - gamma*k'*k(i)/twom;
disp('Clustering ...iterated_genlouvain.m');
%%%%%%Clustering
[SS,QQ,n_it]=iterated_genlouvain(B,limit,0);
max(SS)
toc
[~,r2] = sort(SS);
%         QQ = QQ/twom;
figure;
% imagesc(rCellspks(r2,r2))
AdjComImagesc(rCellspks,SS);
colormap(ColorPN3)
clim([-0.3 0.3])

IndSub=find(SS==3)
CellSub=double(CaData.spks(iscell(IndSub),:)');

[coeff,score,latent,tsquared] = pca(CellSub);
corr(score(:,1:3),fSpeed)



[W,H]=nnmf(double(CaData.spks(iscell,:))',3)
biplot(H','Scores',W);
axis([0 1.1 0 1.1])
xlabel('Column 1')
ylabel('Column 2')

r=corr(W,fSpeed)


% 1 2 3 7 9 13

clear Group
Group(1).Indices=[1:length(tempI)];
MarkPoints3D_GPLmaker(Pos3Dneed, yaml, true, SpiralSizeUM, SpiralRevolutions, SaveName,Group)
MarkPoints3D_XMLmaker_Points(Pos3Dneed,yaml,true, Repetition, SpiralSizeUM, SpiralRevolutions,UncagingLaserPower, SavePathAllPoint)


% MarkPoints3D_XMLmaker_Group(Group,Repetition, SpiralSizeUM, UncagingLaserPower, SavePath)


%  2 4 7 13 14 15 18 23 27

% tempI=[2 3 4 5 9 11 14 16 21 23]
tempI=[2 3 5 9 14 16 21 23]


