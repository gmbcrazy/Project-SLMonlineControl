% Suite2pToSLMResponse(MatFile,BreakPoint,PSTHparam)

MatFile='E:\LuSLMOnlineTest\SL0777-Ai203\11142024\Data\TSeries-11142024-0927-018.mat';
load(MatFile)
% frameRepetion=PreMarkPointRepetition+PostMarkPointRepetition; %%Total repepitions of Z series in T series setting.
nPlane=3;
Ziteration=11;      % NumTrialTseries-1 trials repeated in Single Tseries; As the 1st Zseries is not synchornized with MP. So Actual trials Need to substract 1;
ZRepetition=30;      % number of repetition in each Zseries
frameRepetion=(Ziteration-1)*(ZRepetition-1)+ZRepetition;
TZseriesBin=repmat(30,1,Ziteration);
TZseriesBin=[20 60 40 30 60 20];


TZseriesUpdate=[TZseriesBin(1) TZseriesBin(2:end)-1];
frameBeforeSLM=cumsum(TZseriesUpdate(1:end-1))
frameRepetion=sum(TZseriesUpdate);

PostSLMFrameN=8;
PreSLMFrameN=15;
PointList=[1 1 1 1 2];
CellTarget=[21 37 54];
PointCellList=CellTarget(PointList);

PreSLM=[];
PostSLM=[];
PreSLMMat=[];
PostSLMMat=[];


NeuroData=DeltaF;

for iPoint=1:length(frameBeforeSLM)
    PostStimI=frameBeforeSLM(iPoint)+[1:PostSLMFrameN];
    PreStimI=frameBeforeSLM(iPoint)-PreSLMFrameN+1:frameBeforeSLM(iPoint);

    PreSLM(:,iPoint)=NeuroData(PointCellList(iPoint),PreStimI)';
    PostSLM(:,iPoint)=NeuroData(PointCellList(iPoint),PostStimI)';

    PreSLMMat(:,:,iPoint)=NeuroData(:,PreStimI)';
    PostSLMMat(:,:,iPoint)=NeuroData(:,PostStimI)';


    pRank(iPoint)=ranksum(PreSLM(:,iPoint),PostSLM(:,iPoint));
    [~,pT(iPoint)]=ttest2(PreSLM(:,iPoint),PostSLM(:,iPoint));
end


plot([PreSLM;PostSLM])

figure;
plot(mean([PreSLM;PostSLM],2))


% PSTHMat=cat(1,PreSLMMat,PostSLMMat);
% PSTHaveMat=mean(PSTHMat,3);
% % imagesc(PSTHaveMat')
% 
% figure;
% plot(PSTHaveMat)
% 
% 
% imagesc



