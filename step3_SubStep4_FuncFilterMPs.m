%% Functional Filter of ROIs for SLM targets 
% Including top cells highly correlated associated with speed.
XTimesStdTh = 3;
MinInterVal = 20;
maxLag = 10;
Continueloop=1;


deltaFoF=F2deltaFoF(double(CaData.F),double(CaData.Fneu),double(confSet.fs));
NeuroData=AmpNormalizeRow(deltaFoF',[0 100])';
L=min(size(fSpeed,1),size(NeuroData,1))
NeuroData=NeuroData(1:L,:);
NeuroDataCell=NeuroData(:,iscell);

TopCellNRaw=TopCellN;

for iPlane=1:numPlanes
    I1=find(CaData.CellPlaneID==iPlane);
    [rSpeed(I1),pSpeed(I1)]=corr(NeuroData(:,iscell(I1)),fSpeed(:,iPlane),'type','Spearman','rows','pairwise');
end
clear rStim c
for iCell = 1:size(iscell, 1)
    iCell;
    iPlane=CaData.CellPlaneID(iCell);
    [c(:,iCell), lags] = xcorr(NeuroData(:,iscell(iCell)), fStim(:,iPlane), maxLag, 'coeff');
    PostI = find(lags >= 0);
    [~, i1] = max(abs(c(PostI,iCell)));
    rStim(iCell, 1) = c(PostI(i1),iCell);
end

rSpeed=rSpeed(:);
rStim=rStim(:);

[~,rIStim]=sort(rStim,'descend');  
[~,~,rankStim]=intersect(1:numPoint,rIStim);


[~,rISpeed]=sort(rSpeed,'descend'); 
[~,~,rankSpeed]=intersect(1:numPoint,rISpeed);

%%Cell locates close to edge of the view, were not considered as SML targets.

while Continueloop
    rCenterIStim=intersect(CenterCloseI,rIStim(1:TopCellN));
    rCenterISpeed=intersect(CenterCloseI,rISpeed(1:TopCellN));

    AllFunctionI=union(rCenterIStim,rCenterISpeed);
    ColFunctionI=intersect(rCenterIStim,rCenterISpeed);
    NonFunctionI=setdiff(CenterCloseI,AllFunctionI);

    rCenterIStim=setdiff(rCenterIStim,ColFunctionI);
    rCenterISpeed=setdiff(rCenterISpeed,ColFunctionI);

    if length(rCenterIStim)<TopCellNRaw||length(rCenterISpeed)<TopCellNRaw
        TopCellN=TopCellN+1;
        Continueloop=1;
    else
        Continueloop=0;
    end
end

if length(NonFunctionI)<TopCellN
   NonFunctionI=NonFunctionI;
else
   NonFunctionI=NonFunctionI(randperm(length(NonFunctionI),TopCellN));
end


% plot(lags,c(:,rCenterIStim))
% figure;
% plot(NeuroDataCell(:,rCenterIStim));
% hold on;
% plot(fStim(:,1),'color',[0.4 0.4 0.4])


IncludeCellFunFilter=[rCenterISpeed(:);rCenterIStim(:);NonFunctionI(:)];
IncludeCellFunID=[repmat(1,length(rCenterISpeed),1);repmat(2,length(rCenterIStim),1);repmat(3,length(NonFunctionI),1)];
IncludeCellFunSore=[rSpeed(IncludeCellFunFilter) rStim(IncludeCellFunFilter)];

IncludeInfo=[IncludeCellFunFilter IncludeCellFunID IncludeCellFunSore];
[~,sortI]=sort(IncludeCellFunFilter);
IncludeInfo=IncludeInfo(sortI,:);

IncludeCellFunFilter=IncludeInfo(:,1);
FunScore=IncludeInfo(:,2:end);
SavePathStimSpeed=[SavePath 'Top' num2str(TopCellN) 'SpeedStimEdgeExc\']
mkdir(SavePathStimSpeed)


IncludePathOri=[SavePathStimSpeed '\AllIncludedOrigin\'];
IncludePath=[SavePathStimSpeed '\AllIncluded\'];
mkdir(IncludePathOri)
mkdir(IncludePath)
XYZtoMarkPoint(IncludePathOri,Pos3D,IncludeCellFunFilter,yaml,confSet,CaData.statCell,FunScore);