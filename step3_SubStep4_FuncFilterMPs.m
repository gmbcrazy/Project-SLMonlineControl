%% Functional Filter of ROIs for SLM targets 
% Including top cells highly correlated associated with speed.
XTimesStdTh = 3;
MinInterVal = 20;
maxLag = 10;
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

IncludeCellI=union(rCenterISpeed(1:TopCellN),rCenterIStim(1:TopCellN));
rScore=double([rSpeed(:) rStim(:)]);
rScore=rScore(CenterCloseI,:);
rRank=[rankSpeed(:) rankStim(:)]
rScoreInclude=rScore(IncludeCellI,:);
rRankInclude=rRank(IncludeCellI,:);
IncludeCellFunFilter=CenterCloseI(IncludeCellI);
SavePathStimSpeed=[SavePath 'Top' num2str(TopCellN) 'SpeedStimEdgeExc\']
mkdir(SavePathStimSpeed)


IncludePathOri=[SavePathStimSpeed '\AllIncludedOrigin\'];
IncludePath=[SavePathStimSpeed '\AllIncluded\'];
mkdir(IncludePathOri)
mkdir(IncludePath)
XYZtoMarkPoint(IncludePathOri,Pos3D,IncludeCellFunFilter,yaml,confSet,CaData.statCell);