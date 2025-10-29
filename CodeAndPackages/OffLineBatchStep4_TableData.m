clear all


BatchSavePath='D:\Project1-LocalProcessing\Step1\';
load([BatchSavePath '07-Oct-2025FOV.mat'])
Suite2pDataKeywords='awakeRefSpon';

DataSavePath='\\nimhlabstore1.nimh.nih.gov\UFNC\FNC2\Zhang\Projects\Project-LocalProcessing\Step4\';
mkdir(DataSavePath);
DataSavePath=[DataSavePath Suite2pDataKeywords '\'];
mkdir(DataSavePath);
% SaveFunCon=[DataSavePath 'FunCon\'];
% mkdir(SaveFunCon)

% VolEventLabel={}

PSTHparam.PreSLMCal = 10; 
PSTHparam.PostSLMCal = 15;
PSTHparam.pTh = 0.05; 
PSTHparam.TestMethod = 'ranksum';
PSTHparam.MPFrameJump = 2;
PSTHparam.TestStepFrame = 3;    %%post-slm frames for Test whether SLM works
PSTHparam.PreTestFrame = 10; 
PSTHparam.iData = 1;    %%post-slm frames for Test whether SLM works
umPerlPixel=700/512;
offTargetum=15;
PSTHparam.OffTargetPixel=ceil(offTargetum/umPerlPixel)

 


ResponseMap=slanCM('seismic',64);

TimBinFrame = -PSTHparam.PreSLMCal:PSTHparam.PostSLMCal-1;
IndexFOVNeed=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20];



for iData=1:2
    PSTHparam.iData = iData;    %%post-slm frames for Test whether SLM works
    [Output(iData),NeuroTrace,BehTrace]=OfflineSLM_ExtractFOVs(FOVUpdate(IndexFOVNeed), Suite2pDataKeywords,suite2pFOVPath(IndexFOVNeed),PSTHparam);
end


CellN=size(Output(1).NeuroPos3DMeta,1);
CellNNonTarget=size(Output(1).NeuroPos3DMeta,1)-size(Output(1).GroupTargetCellAll,1);


% SaveFunCon=[DataSavePath 'GroupSLM' num2str(length(IndexFOVNeed)) 'Sessions\'];
% mkdir(SaveFunCon)

% SaveFunCon=[DataSavePath 'GroupSLM' num2str(length(IndexFOVNeed)) 'Sessions\'];
SaveFunFOV=[DataSavePath 'GroupSLM' num2str(length(IndexFOVNeed)) 'Sessions\'];
mkdir(SaveFunFOV)


save([SaveFunFOV 'FOVoutputs.mat'],'-v7.3');
load([SaveFunFOV 'FOVoutputs.mat'])


PostPreDiffSpeedTh=[1 2 10000];
iSpeedTh=2;
ProcessPar.PostPreDiffSpeedTh=PostPreDiffSpeedTh(iSpeedTh);    %%noticed that this speed threshold is not applied when use OfflineSLM_FOVmeta2NeuroDeltaTrial.
ProcessPar.PlotFigure=0;
ProcessPar.GroupLabel={'L','S','N'};
ProcessPar.GroupList=[1 2 3];
ProcessPar.GroupColor=[255 51 153;91 20 212;121 247 111]/255;
ProcessPar.PowerZeroColor=[0.5 0.5 0.5];
ProcessPar.PowerZero=[0 1];
ProcessPar.PowerZeroLabel={'SLM','FakeSLM'};
ProcessPar.VolOut=[0 1];
ProcessPar.VolOutLabel={'NoWhisk','Whisk'}
ProcessPar.AwakeState=[1];
ProcessPar.AwakeStateLabel={'Awake'};
%%Making table for mixed linear models with cell and trial as sample
NDataName={'deltaF','spks'};

for iData=1:2

    PSTHparam.iData=iData;

PSTHparam.TestWinNum=5;
tblResponse=[];
for iWin=1:PSTHparam.TestWinNum
    PSTHparam.PostWinN=iWin;
    [tbl,InfoTable]=OfflineSLM_FOVmeta2NeuroDeltaTrial(Output(iData),ProcessPar,PSTHparam);

    % SLMGroupTableTrial=[tbl InfoTable];
    NeedInfo={'Group','VolOut','PowerZero','AwakeState'};

    SLMGroupTableTrial = [tbl InfoTable(:,NeedInfo)];
    SLMGroupTableTrial(isnan(SLMGroupTableTrial.Group),:)=[];
    SLMGroupTableTrial=removevars(SLMGroupTableTrial,{'PointSpeedR','PointStimR'});
    
    SLMGroupTableTrial.Properties.VariableNames{'VolOut'} = 'Whisk';
    
    SLMGroupTrialID=SLMGroupTableTrial.TrialID;
    EndSessInd=[find(diff(SLMGroupTrialID)<0);length(SLMGroupTrialID)];
    StartSessInd=[1;find(diff(SLMGroupTrialID)<0)+1];
    Time=[];
    for iT=1:length(EndSessInd)
        TempSeq=StartSessInd(iT):EndSessInd(iT);
        TempSeq=(TempSeq-min(TempSeq))/(max(TempSeq)-min(TempSeq));
        Time=[Time;TempSeq(:)];
    end
    SLMGroupTableTrial.Time=Time;
    writetable(SLMGroupTableTrial,[SaveFunFOV NDataName{iData} 'SLMGroupResTrialDynWin' num2str(iWin) '.csv']);
    tblResponse(:,iWin)=tbl.Response;
end

end

PSTHparam.TestWinNum=1;

for iData=1:2

    PSTHparam.iData=iData;

tblResponse=[];
for iWin=1:PSTHparam.TestWinNum
    PSTHparam.PostWinN=iWin;
    [tbl,InfoTable]=OfflineSLM_FOVmeta2NeuroDeltaTrial(Output(iData),ProcessPar,PSTHparam);
    % SLMGroupTableTrial=[tbl InfoTable];
    NeedInfo={'UncagingLaserPower','Point','PointTargetCell','PointTargetCellGroup'};

    SLMPointTableTrial = [tbl InfoTable(:,NeedInfo)];
    SLMPointTableTrial(isnan(SLMPointTableTrial.Point),:)=[];
    % SLMPointTableTrial=removevars(SLMPointTableTrial,{'PointSpeedR','PointStimR'});
    
    % SLMPointTableTrial.Properties.VariableNames{'VolOut'} = 'Whisk';
    
    SLMPointTrialID=SLMPointTableTrial.TrialID;
    EndSessInd=[find(diff(SLMPointTrialID)<0);length(SLMPointTrialID)];
    StartSessInd=[1;find(diff(SLMPointTrialID)<0)+1];
    Time=[];
    for iT=1:length(EndSessInd)
        TempSeq=StartSessInd(iT):EndSessInd(iT);
        TempSeq=(TempSeq-min(TempSeq))/(max(TempSeq)-min(TempSeq));
        Time=[Time;TempSeq(:)];
    end
    SLMPointTableTrial.Time=Time;
    writetable(SLMPointTableTrial,[SaveFunFOV NDataName{iData} 'SLMPointResTrialDynWin' num2str(iWin) '.csv']);
    tblResponse(:,iWin)=tbl.Response;
end

end



SLMGroupTableTrial=readtable([SaveFunFOV 'deltaFSLMGroupResTrialDynWin1.csv']);



DimentionReductionMethod='NonNegMatFac';
SaveFunFOVDim=[SaveFunFOV DimentionReductionMethod '\'];
mkdir(SaveFunFOVDim)


SponInd=[1:8200];
ExcludeStimInd=[1:2000 2800:4100];  
ExcludeStimInd=union(ExcludeStimInd,ExcludeStimInd+4100);

% AddStimNoise=random('norm',0,0.1,length(SponInd),1);

SaveFunSub=[SaveFunFOVDim 'FOV\'];
mkdir(SaveFunSub)
% clear FOVres

SaveFunFOVDim
iData=1;
pcaN=40;
ScoreLabel={'SpeedScore','WhiskScore'}


% for iFOV=3:3

for iFOV=1:length(IndexFOVNeed)
    close all
tblFOV=SLMGroupTableTrial(SLMGroupTableTrial.Session==iFOV,:);
TrialID=unique(tblFOV.TrialID);
tbltemp=SLMGroupTableTrial(SLMGroupTableTrial.Session==iFOV&SLMGroupTableTrial.TrialID==TrialID(1),:);
% [W,H] = nnmf(NeuroTrace{iFOV}{iData}(tbltemp.TargetCell==0,SponInd),pcaN);
[W,H] = nnmf(NeuroTrace{iFOV}{iData}(:,SponInd),pcaN,'algorithm','als','replicates',20);

I1=find(tblFOV.TrialID==TrialID(1));
NonTarget=find(tblFOV.TargetCell(I1)==0);
Target=find(tblFOV.TargetCell(I1)==1);

% Wsub=W(NonTarget,:);
% norms = sqrt( sum( Wsub.^2 , 1 ) );  
% WNontarget = bsxfun(@rdivide, Wsub, norms);  
% Asub=NeuroTrace{iFOV}{iData}(NonTarget,SponInd);

WNontarget=W(NonTarget,:);
% end
pinvW = pinv(W);
subNonTargetData=NeuroTrace{iFOV}{iData}(:,SponInd);
% subNonTargetData(FOVres(iFOV).Target,:)=0;
subNonTargetData(Target,:)=0;
HNontarget = (pinvW * subNonTargetData); 
FOVres(iFOV).W=W;
FOVres(iFOV).H=H;
FOVres(iFOV).WNontarget=WNontarget;
FOVres(iFOV).HNontarget=HNontarget;
FOVres(iFOV).Nontarget=NonTarget;
FOVres(iFOV).Target=Target;




SubData=NeuroTrace{iFOV}{iData}(NonTarget,SponInd);
[W_sub,H_sub] = nnmf(SubData,pcaN,'algorithm','als','replicates',20);
FOVres(iFOV).W_sub=W_sub;
FOVres(iFOV).H_sub=H_sub;

% HNontarget=WNontarget'*NeuroTrace{iFOV}{iData}(NonTarget,SponInd);

% for t=1:length(SponInd)
%     HNontarget(:,t) = lsqnonneg(WNontarget, SubData(:,t));



% nComp = size(WNontarget,2);
% nCols = size(Asub,2);
% Hnorm = zeros(nComp,nCols);
% for j = 1:nCols
%   Hnorm(:,j) = lsqnonneg(WNontarget, Asub(:,j) );
% end
% HNontarget = bsxfun(@times, Hnorm, norms(:));
% HNontarget=Hnorm;
% 
% H=H';
% HNontarget=HNontarget';
ProcessPar.GroupColor=[255 51 153;91 20 212;121 247 111]/255;

% close all
% for PCI=1:10
%     figure;
%     [~,rankI]=sort((W(:,PCI)),'descend');
%     imagesc(NeuroTrace{iFOV}{iData}(rankI,SponInd));
%     % subplot(1,2,1)
%     % plot(AmpNormalize(coeff(SponInd,PCI),[1 99]),'k');
%     % hold on;plot(AmpNormalize(nanmean(BehTrace(iFOV).Speed(SponInd,:),2)),'g')
%     % % plot(AmpNormalize(coeff(SponInd,1),[1 99]),'r');
%     % hold on;plot(AmpNormalize(nanmean(BehTrace(iFOV).Stim(SponInd,:),2)),'r')
% end
% 
HH=H';

rSpeed=corr(HH(ExcludeStimInd,1:pcaN),nanmean(BehTrace(iFOV).Speed(ExcludeStimInd,:),2),'rows','complete');


% rSpeed=corr(H(:,1:pcaN),nanmean(BehTrace(iFOV).Speed(:,:),2),'rows','complete');

% rSpeedNon=corr(H(ExcludeStimInd,1:pcaN),nanmean(BehTrace(iFOV).Speed(ExcludeStimInd,:),2),'rows','complete');
maxLag=10;
clear rStim c rStim_sub c_sub


for PCI = 1:pcaN
    [c(:,PCI), lags] = xcorr(AmpNormalize(HH(SponInd,PCI),[0 100]), AmpNormalize(nanmean(BehTrace(iFOV).Stim(SponInd,:),2),[0 100]), maxLag, 'coeff');
   
    PostI = find(lags >= 0&lags <= 10);
    PreI=find(lags<0);
    % rStim(PCI, 1)=max(c(PostI,PCI));
    % [~, i1] = max(c(PostI,PCI));
    
    % [~, i1] = max(abs(c(PostI,PCI)-mean(c(PreI,PCI))));

    %%Current used
    [~, i1] = max(abs(c(PostI,PCI))-mean(c(PreI,PCI)));
    rStim(PCI, 1) = c(PostI(i1),PCI)-mean(AmpNormalize(HH(SponInd,PCI)));
    %%Current used


end
[~,pcaIspeed]=max(rSpeed)
[~,pcaIstim]=max(rStim)

FOVres(iFOV).SpeedScore=rSpeed([pcaIspeed pcaIstim]);
FOVres(iFOV).StimScore=rStim([pcaIspeed pcaIstim]);


% for PCI = 1:pcaN
%     figure
%     plot(AmpNormalize(HH_sub(SponInd,PCI)));
%     hold on;
%     plot(nanmean(BehTrace(iFOV).Stim(SponInd,:),2)>1,'color',ProcessPar.GroupColor(2,:))
% end


HH_sub=H_sub';
SpeedSpon=nanmean(BehTrace(iFOV).Speed(ExcludeStimInd,:),2);
rSpeed_sub=corr(HH_sub(ExcludeStimInd,1:pcaN),SpeedSpon,'rows','complete');

for PCI = 1:pcaN
    [c_sub(:,PCI), lags] = xcorr(AmpNormalize(HH_sub(SponInd,PCI),[0 100]), AmpNormalize(nanmean(BehTrace(iFOV).Stim(SponInd,:),2),[0 100]), maxLag, 'coeff');
    PostI = find(lags >= 0&lags <= 10);
    PreI=find(lags<0);
  %%Current used
    [~, i1] = max(abs(c_sub(PostI,PCI))-mean(c_sub(PreI,PCI)));
    rStim_sub(PCI, 1) = c_sub(PostI(i1),PCI)-mean(AmpNormalize(HH_sub(SponInd,PCI)));
    %%Current used


end

[~,pcaIspeed_sub]=max(rSpeed_sub)
[~,pcaIstim_sub]=max(rStim_sub)

FOVres(iFOV).SpeedScore_sub=rSpeed_sub([pcaIspeed_sub pcaIstim_sub]);
FOVres(iFOV).StimScore_sub=rStim_sub([pcaIspeed_sub pcaIstim_sub]);


[FOVres(iFOV).SpeedReg.B,FOVres(iFOV).SpeedReg.BINT,FOVres(iFOV).SpeedReg.R,...
    FOVres(iFOV).SpeedReg.RINT,FOVres(iFOV).SpeedReg.STATS]=regress(H(pcaIspeed,ExcludeStimInd)',[ones(length(SpeedSpon),1) SpeedSpon]);

[FOVres(iFOV).SpeedRegNontarget.B,FOVres(iFOV).SpeedRegNontarget.BINT,FOVres(iFOV).SpeedRegNontarget.R,...
    FOVres(iFOV).SpeedRegNontarget.RINT,FOVres(iFOV).SpeedRegNontarget.STATS]=regress(HNontarget(pcaIspeed,ExcludeStimInd)',[ones(length(SpeedSpon),1) SpeedSpon]);

[FOVres(iFOV).SpeedReg_sub.B,FOVres(iFOV).SpeedReg_sub.BINT,FOVres(iFOV).SpeedReg_sub.R,...
    FOVres(iFOV).SpeedReg_sub.RINT,FOVres(iFOV).SpeedReg_sub.STATS]=regress(H_sub(pcaIspeed_sub,ExcludeStimInd)',[ones(length(SpeedSpon),1) SpeedSpon]);

% x=[min(data1(:)) max(data1(:))]';
% y=[[1 1]' x]*B;
% close all
% plot(lags,c_sub(:,12))

pause(3)

FOVres(iFOV).ComI=[pcaIspeed pcaIstim];
FOVres(iFOV).ComI_sub=[pcaIspeed_sub pcaIstim_sub];


WSub=W(:,[pcaIspeed pcaIstim]);
WNontargetSub=WNontarget(:,[pcaIspeed pcaIstim]);

PP=[0.1 0.03 0.05 0.1];
      % P.xLeft=PP(1);
      % P.xRight=PP(2);
      % P.yTop=PP(3);
      % P.yBottom=PP(4);


close all
% for pcaIspeed=1:11
figure;
PosMatX=[0.85;0.85;0.85];
PosMatY=[0.25;0.25;0.25];
t=subplotPosLu(PosMatX,PosMatY,1,1,PP);
plot(H(pcaIspeed,SponInd),'Color',[0.1 0.1 0.1]);hold on;
plot(HNontarget(pcaIspeed,SponInd),'Color',[0.7 0.7 0.7]);
text(7000, max(H(pcaIspeed,SponInd)),num2str(FOVres(iFOV).SpeedScore(1)),'Color',ProcessPar.GroupColor(1,:),'horizontalalignment','right')
text(8200, max(H(pcaIspeed,SponInd)),num2str(FOVres(iFOV).StimScore(1)),'Color',ProcessPar.GroupColor(2,:),'horizontalalignment','right')

ylabel(['AllNeuroCom.' num2str(pcaIspeed)]);
set(gca,'xlim',[0 max(SponInd)],'box','off','xticklabel',[]);

t=subplotPosLu(PosMatX,PosMatY,2,1,PP);
plot(H_sub(pcaIspeed_sub,SponInd),'Color',[0.1 0.1 0.1]);
set(gca,'xlim',[0 max(SponInd)],'box','off','xticklabel',[]);
ylabel(['NTNeuroCom.' num2str(pcaIspeed_sub)]);
set(gca,'xlim',[0 max(SponInd)],'box','off','xticklabel',[]);
text(7000, max(H_sub(pcaIspeed_sub,SponInd)),num2str(FOVres(iFOV).SpeedScore_sub(1)),'Color',ProcessPar.GroupColor(1,:),'horizontalalignment','right')
text(8200, max(H_sub(pcaIspeed_sub,SponInd)),num2str(FOVres(iFOV).StimScore_sub(1)),'Color',ProcessPar.GroupColor(2,:),'horizontalalignment','right')

subplotPosLu(PosMatX,PosMatY,3,1,PP)
plot(nanmean(BehTrace(iFOV).Speed(SponInd,:),2),'color',ProcessPar.GroupColor(1,:))
ylabel(['Speed'])
set(gca,'xlim',[0 max(SponInd)],'box','off');
xlabel('Frames')
% plot(AmpNormalize(coeff(SponInd,1),[1 99]),'r');
% hold on;plot(AmpNormalize(nanmean(BehTrace(iFOV).Stim(SponInd,:),2)),'r')
% end
    papersizePX=[0 0 15 10];
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
    print(gcf, [SaveFunSub  NDataName{iData} 'FOV' num2str(iFOV) 'NeuroComSpeed.svg'], '-dsvg', '-painters');
    print(gcf, [SaveFunSub  NDataName{iData} 'FOV' num2str(iFOV) 'NeuroComSpeed.tif'], '-dtiffn', '-painters');


figure;
t=subplotPosLu(PosMatX,PosMatY,1,1,PP)
plot(H(pcaIstim,SponInd),'Color',[0.1 0.1 0.1]);hold on;
plot(HNontarget(pcaIstim,SponInd),'Color',[0.7 0.7 0.7]);
set(gca,'xlim',[0 max(SponInd)],'box','off','xticklabel',[]);
ylabel(['NeuroCom.' num2str(pcaIstim)])
text(7000, max(H(pcaIstim,SponInd)),num2str(FOVres(iFOV).SpeedScore(2)),'Color',ProcessPar.GroupColor(1,:),'horizontalalignment','right')
text(8200, max(H(pcaIstim,SponInd)),num2str(FOVres(iFOV).StimScore(2)),'Color',ProcessPar.GroupColor(2,:),'horizontalalignment','right')


t=subplotPosLu(PosMatX,PosMatY,2,1,PP)
plot(H_sub(pcaIstim_sub,SponInd),'Color',[0.1 0.1 0.1]);
ylabel(['NTNeuroCom.' num2str(pcaIstim_sub)]);
set(gca,'xlim',[0 max(SponInd)],'box','off','xticklabel',[]);
text(7000, max(H_sub(pcaIstim_sub,SponInd)),num2str(FOVres(iFOV).SpeedScore_sub(2)),'Color',ProcessPar.GroupColor(1,:),'horizontalalignment','right')
text(8200, max(H_sub(pcaIstim_sub,SponInd)),num2str(FOVres(iFOV).StimScore_sub(2)),'Color',ProcessPar.GroupColor(2,:),'horizontalalignment','right')

subplotPosLu(PosMatX,PosMatY,3,1,PP)
plot(nanmean(BehTrace(iFOV).Stim(SponInd,:),2)>1,'color',ProcessPar.GroupColor(2,:))
ylabel(['Whisker'])
set(gca,'xlim',[0 max(SponInd)],'box','off');
xlabel('Frames')

% plot(AmpNormalize(coeff(SponInd,1),[1 99]),'r');
    papersizePX=[0 0 15 10];
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
    print(gcf, [SaveFunSub NDataName{iData} 'FOV' num2str(iFOV) 'NeuroComWhisker.svg'], '-dsvg', '-painters');
    print(gcf, [SaveFunSub NDataName{iData} 'FOV' num2str(iFOV) 'NeuroComWhisker.tif'], '-dtiffn', '-painters');


end

TBLall=[];
TBLallNontarget=[];
TBLall_sub=[];
    
for iFOV=1:length(IndexFOVNeed)
close all
tblFOV=SLMGroupTableTrial(SLMGroupTableTrial.Session==iFOV,:);
TrialID=unique(tblFOV.TrialID);


DataVecInfo=[];
DataVec=[];
ScoreNontarget=[];
Score=[];

Score_sub=[];

% DataVecTBL=[];
clear DataVecTBL
tempScore=[];
for iTrial=1:length(TrialID)
    DataVec(:,iTrial)= tblFOV.Response(tblFOV.TrialID==TrialID(iTrial),:);
    I1=find(tblFOV.TrialID==TrialID(iTrial));

    % I1=find(tblFOV.TrialID==TrialID(1));
    NonTarget=find(tblFOV.TargetCell(I1)==0);
    Target=find(tblFOV.TargetCell(I1)==1);

    temp=tblFOV(tblFOV.TrialID==TrialID(iTrial),:);
    % DataVecInfo(iTrial,:)=[tblFOV.Group(I1(1)) tblFOV.Whisk(I1(1)) tblFOV.PowerZero(I1(1)) tblFOV.Time(I1(1))  tblFOV.Speed(I1(1))];
    % DataVec(tblFOV.TargetCell(I1)==1,iTrial)=0;
    DataVecTBL(iTrial,:)=tblFOV(I1(1),:);
    DataVecTBL.SpeedR(iTrial)=mean(tblFOV.SpeedR(I1(NonTarget),:));
    DataVecTBL.StimR(iTrial)=mean(tblFOV.StimR(I1(NonTarget),:));
    % DataVecTBL(iTrial,:).TargetSpeedR=mean(tblFOV.SpeedR(I1(Target),:));
    % DataVecTBL(iTrial,:).TargetStimR=mean(tblFOV.StimR(I1(Target),:));
    if temp.Group(1)~=3
        tempScore(iTrial,:)=[FOVres(iFOV).SpeedScore(temp.Group(1)) FOVres(iFOV).StimScore(temp.Group(1))];
        % DataVecTBL.TargetStimScore(iTrial)=FOVres.StimScore(temp.Group(1));     
    else
        tempScore(iTrial,:)=[NaN NaN];
    end


    Score(:,iTrial)=lsqnonneg(FOVres(iFOV).W(:,FOVres(iFOV).ComI), DataVec(:,iTrial));
    ScoreNontarget(:,iTrial)=lsqnonneg(FOVres(iFOV).WNontarget(:,FOVres(iFOV).ComI), DataVec(NonTarget,iTrial));
    Score_sub(:,iTrial)=lsqnonneg(FOVres(iFOV).W_sub(:,FOVres(iFOV).ComI_sub), DataVec(NonTarget,iTrial));

end

DataVecTBL=[DataVecTBL array2table(tempScore,"VariableNames",{'TargetSpeedScore','TargetWhiskerScore'})]
Score=Score';
ScoreNontarget=ScoreNontarget';
Score_sub=Score_sub';


DataVecRaw=DataVec;
DataVecNontarget=DataVec(NonTarget,:);

% Score=DataVec'*WSub;
% ScoreNontarget=DataVecNontarget'*WNontargetSub;

FOVres(iFOV).DataVec=DataVec;
FOVres(iFOV).InfoTBL=DataVecTBL;
% FOVres(iFOV).NonTarget=NonTarget;
FOVres(iFOV).DataVecNonTarget=DataVecNontarget;


    pinvW = pinv(FOVres(iFOV).W(:,FOVres(iFOV).ComI));
    Score = (pinvW * FOVres(iFOV).DataVec)'; 
    Score(:,1) = Score(:,1) - [ones(length(Score(:,1)),1) DataVecTBL.Speed]*FOVres(iFOV).SpeedReg.B;


    % pinvW = pinv(FOVres(iFOV).W(:,FOVres(iFOV).ComI));
    subNonTargetData=FOVres(iFOV).DataVec;
    % subNonTargetData(FOVres(iFOV).Target,:)=0;
    subNonTargetData(Target,:)=0;
    % ScoreNontarget = (pinvW * subNonTargetData)'; 
    ScoreNontarget(:,1) = ScoreNontarget(:,1) - [ones(length(ScoreNontarget(:,1)),1) DataVecTBL.Speed]*FOVres(iFOV).SpeedRegNontarget.B;

    pinvW_sub = pinv(FOVres(iFOV).W_sub(:,FOVres(iFOV).ComI_sub));
    % Score_sub = (pinvW_sub * FOVres(iFOV).DataVecNonTarget)'; 
    Score_sub(:,1) = Score_sub(:,1) - [ones(length(Score_sub(:,1)),1) DataVecTBL.Speed]*FOVres(iFOV).SpeedReg_sub.B;

maxScore=max(Score);
% maxScore=max(ScoreNontarget)

% SpeedTh=0.5;
clear ScoreGroup SpeedGroup ;
% maxScore
for iWhisk=1:2
for iGroup=1:3
    I1=find(DataVecTBL.Group==iGroup&DataVecTBL.Whisk==ProcessPar.VolOut(iWhisk)&DataVecTBL.PowerZero==0&abs(DataVecTBL.Speed)<1000);
    ScoreGroup(iGroup,iWhisk).Score=Score(I1,:);
    ScoreGroup(iGroup,iWhisk).ScoreNonTarget=ScoreNontarget(I1,:);
    ScoreGroup(iGroup,iWhisk).Time=DataVecTBL.Time(I1);
    SpeedGroup{iGroup,iWhisk}=DataVecTBL.Speed(I1);
    ScoreGroup(iGroup,iWhisk).Score_sub=Score_sub(I1,:);
end

% figure;hold on;



figure;
for iScore=1:2
TempScore={};

    for iGroup=1:3
        % plot(ScoreGroup(iGroup).Time,ScoreGroup(iGroup).Score(:,1),'Color',ProcessPar.GroupColor(iGroup,:));
        NeuroScore{iGroup,iScore}=ScoreGroup(iGroup,iWhisk).Score(iScore,:);
        NeuroScoreNontarget{iGroup,iScore}=ScoreGroup(iGroup,iWhisk).ScoreNonTarget(iScore,:);
        NeuroScore_sub{iGroup,iScore}=ScoreGroup(iGroup,iWhisk).Score_sub(iScore,:);


    end


   GroupPair.CorrName='fdr';
   GroupPair.Q=0.1;
   GroupPair.Pair=[];
   GroupPair.SignY=maxScore(iScore);
   GroupPair.Plot=1;
   GroupPair.Std=0;      %%%%%%%%%using standard deviation as errorbar
   GroupPair.SamplePlot=1; %%%%%%%%%Plot Individual Sample Point
   GroupPair.SamplePairedPlot=0; %%%%%%%%%Dash line for paired comparison sample
   GroupPair.LimY=[-0.005 GroupPair.SignY*1.6];
   GroupPair.Marker={'o'};   % Datatype=varargin{1};  %%%%%%Datatype=0, non-paired data; 1 for paired data
   % SavePath=varargin{2};  %%%%%%Save the stats to txt file
   % GroupPair=varargin{3};
   % RGroupID=varargin{4};


   GroupPairSpeed=GroupPair;
   GroupPairSpeed.SignY=3;
   GroupPairSpeed.LimY=[-0.2 GroupPairSpeed.SignY*1.2];

% stats=ErrorBarPlotLU(1:3,NeuroScore,[],ProcessPar.GroupColor,2,0,[],GroupPair)
   % SavePath=varargin{4};
   % GroupPair=varargin{5};
   %
    subplotLU(2,4,iScore,1)
    % if size(NeuroScore,2)<3
    %    continue;
    % end
    stats=ErrorBoxPlotLU(1:3,NeuroScore(:,iScore),ProcessPar.GroupColor,[],GroupPair,[1 2 3]);
    set(gca,'ylim',GroupPair.LimY);

    if iScore==2
       xlabel('All')
    end

     ylabel(ScoreLabel{iScore})



    subplotLU(2,4,iScore,2)
    stats=ErrorBoxPlotLU(1:3,NeuroScoreNontarget(:,iScore),ProcessPar.GroupColor,[],GroupPair,[1 2 3])
    set(gca,'ylim',GroupPair.LimY);

    if iScore==2
       xlabel('NonTarget')
    end

    subplotLU(2,4,iScore,3)
    stats=ErrorBoxPlotLU(1:3,NeuroScore_sub(:,iScore),ProcessPar.GroupColor,[],GroupPair,[1 2 3])
    set(gca,'ylim',GroupPair.LimY);

    if iScore==2
       xlabel('NonTarget')
    end


    if iScore==2
    subplotLU(2,4,iScore,4)
    stats=ErrorBoxPlotLU(1:3,SpeedGroup(:,iWhisk),ProcessPar.GroupColor,[],GroupPairSpeed,[1 2 3])
    set(gca,'ylim',GroupPairSpeed.LimY);

    % if iScore==2
       xlabel('Speed')
    % end
    end
    % subplotLU(2,2,1,1)
    % stats=ErrorBoxPlotLU(1:3,NeuroScore(:,iScore),ProcessPar.GroupColor,[],GroupPair,[1 2 3]);
    % subplotLU(2,2,1,2)
    % stats=ErrorBoxPlotLU(1:3,NeuroScoreNontarget(:,iScore),ProcessPar.GroupColor,[],GroupPair,[1 2 3])

end


    papersizePX=[0 0 15 12];
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
    print(gcf, [SaveFunSub ProcessPar.VolOutLabel{iWhisk} 'FOV' num2str(iFOV) 'NeuroCom' '.svg'], '-dsvg', '-painters');
    print(gcf, [SaveFunSub ProcessPar.VolOutLabel{iWhisk} 'FOV' num2str(iFOV) 'NeuroCom' '.tif'], '-dtiffn', '-painters');



end



    ScoreNontarget=array2table(ScoreNontarget,"VariableNames",{'SpeedScore','WhiskerScore'});
    Score=array2table(Score,"VariableNames",{'SpeedScore','WhiskerScore'});
    Score_sub=array2table(Score_sub,"VariableNames",{'SpeedScore','WhiskerScore'});

    TBLallNontarget=[TBLallNontarget;[ScoreNontarget FOVres(iFOV).InfoTBL(:,2:end)]];
    TBLall=[TBLall;[Score FOVres(iFOV).InfoTBL(:,2:end)]];
    TBLall_sub=[TBLall_sub;[Score_sub FOVres(iFOV).InfoTBL(:,2:end)]];





clear DataVecTBL DataVec DataVecRaw DataVecNontarget

end



save([SaveFunSub NDataName{iData} date 'nnmfResult.mat'],'-v7.3');



TBLall=[];
TBLallNontarget=[];
TBLall_sub=[];

for iFOV=1:length(IndexFOVNeed)
    % if iFOV==6
    %    continue;
    % end

    tblFOV=SLMGroupTableTrial(SLMGroupTableTrial.Session==iFOV,:);
    TrialID=unique(tblFOV.TrialID);
    ScoreNontarget=[];
    Score=[];
    Score_sub=[];
    DataVec=[];
    tempScore=[];
    clear DataVecTBL
    for iTrial=1:size(FOVres(iFOV).DataVec,2)
        DataVec(:,iTrial)= tblFOV.Response(tblFOV.TrialID==TrialID(iTrial),:);
        I1=find(tblFOV.TrialID==TrialID(iTrial));
    
        % I1=find(tblFOV.TrialID==TrialID(1));
        NonTarget=find(tblFOV.TargetCell(I1)==0);
        Target=find(tblFOV.TargetCell(I1)==1);
    
        % temp=tblFOV(tblFOV.TrialID==TrialID(iTrial),:);
        % DataVecInfo(iTrial,:)=[tblFOV.Group(I1(1)) tblFOV.Whisk(I1(1)) tblFOV.PowerZero(I1(1)) tblFOV.Time(I1(1))  tblFOV.Speed(I1(1))];
        % DataVec(tblFOV.TargetCell(I1)==1,iTrial)=0;
        DataVecTBL(iTrial,:)=tblFOV(I1(1),:);
        DataVecTBL.SpeedR(iTrial)=mean(tblFOV.SpeedR(I1(NonTarget),:));
        DataVecTBL.StimR(iTrial)=mean(tblFOV.StimR(I1(NonTarget),:));
        % DataVecTBL(iTrial,:).TargetSpeedR=mean(tblFOV.SpeedR(I1(Target),:));
        % DataVecTBL(iTrial,:).TargetStimR=mean(tblFOV.StimR(I1(Target),:));
        if temp.Group(1)~=3
            tempScore(iTrial,:)=[FOVres(iFOV).SpeedScore(temp.Group(1)) FOVres(iFOV).StimScore(temp.Group(1))];
            % DataVecTBL.TargetStimScore(iTrial)=FOVres.StimScore(temp.Group(1));     
        else
            tempScore(iTrial,:)=[NaN NaN];
        end
        % 
        % Keep Score positive
        Score(:,iTrial)=lsqnonneg(FOVres(iFOV).W(:,FOVres(iFOV).ComI), FOVres(iFOV).DataVec(:,iTrial));
        ScoreNontarget(:,iTrial)=lsqnonneg(FOVres(iFOV).WNontarget(:,FOVres(iFOV).ComI), FOVres(iFOV).DataVecNonTarget(:,iTrial));
        Score_sub(:,iTrial)=lsqnonneg(FOVres(iFOV).W_sub(:,FOVres(iFOV).ComI_sub), FOVres(iFOV).DataVecNonTarget(:,iTrial));

        % %%Keep Score positive

    end

    DataVecTBL=[DataVecTBL array2table(tempScore,"VariableNames",{'TargetSpeedScore','TargetWhiskerScore'})];

    Score=Score';
    ScoreNontarget=ScoreNontarget';
    Score_sub=Score_sub';
    % 
    % Score=zscore(Score);
    % ScoreNontarget=zscore(ScoreNontarget);
    % Score_sub=zscore(Score_sub);
    % 
    % pinvW = pinv(FOVres(iFOV).W(:,FOVres(iFOV).ComI));
    % Score = (pinvW * FOVres(iFOV).DataVec)'; 
    Score(:,1) = Score(:,1) - [ones(length(Score(:,1)),1) DataVecTBL.Speed]*FOVres(iFOV).SpeedReg.B;

% x=[min(data1(:)) max(data1(:))]';
% y=[[1 1]' x]*B;
% close all
% plot(lags,c_sub(:,12))
    % subNonTargetData=FOVres(iFOV).DataVec;
    % subNonTargetData(Target,:)=0;
    % ScoreNontarget = (pinvW * subNonTargetData)'; 
    ScoreNontarget(:,1) =  ScoreNontarget(:,1) - [ones(length(Score(:,1)),1) DataVecTBL.Speed]*FOVres(iFOV).SpeedRegNontarget.B;

    % pinvW_sub = pinv(FOVres(iFOV).W_sub(:,FOVres(iFOV).ComI_sub));
    % Score_sub = (pinvW_sub * FOVres(iFOV).DataVecNonTarget)'; 
    % 
    Score_sub(:,1) = Score_sub(:,1) - [ones(length(Score(:,1)),1) DataVecTBL.Speed]*FOVres(iFOV).SpeedReg_sub.B;

    % ScoreNontarget=array2table(ScoreNontarget,"VariableNames",{'SpeedScore','WhiskerScore'});
    % Score=array2table(Score,"VariableNames",{'SpeedScore','WhiskerScore'});

    % ScoreNontarget=zscore(FOVres(iFOV).DataVecNonTarget'*FOVres(iFOV).WNontarget(:,FOVres(iFOV).ComI));
    % ScoreNontarget=array2table(ScoreNontarget',"VariableNames",{'SpeedScore','WhiskerScore'});
    % Score=array2table(Score',"VariableNames",{'SpeedScore','WhiskerScore'});
    % 
    ScoreNontarget=array2table(ScoreNontarget,"VariableNames",{'SpeedScore','WhiskerScore'});
    Score=array2table(Score,"VariableNames",{'SpeedScore','WhiskerScore'});
    Score_sub=array2table(Score_sub,"VariableNames",{'SpeedScore','WhiskerScore'});


    % ScoreNontarget=array2table(zscore(ScoreNontarget),"VariableNames",{'SpeedScore','WhiskerScore'});
    % Score=array2table(zscore(Score),"VariableNames",{'SpeedScore','WhiskerScore'});
    % Score_sub=array2table(zscore(Score_sub),"VariableNames",{'SpeedScore','WhiskerScore'});

    % ScoreNontarget=array2table(ScoreNontarget,"VariableNames",{'SpeedScore','WhiskerScore'});
    % Score=array2table(Score,"VariableNames",{'SpeedScore','WhiskerScore'});

    TBLallNontarget=[TBLallNontarget;[ScoreNontarget FOVres(iFOV).InfoTBL(:,2:end)]];
    TBLall=[TBLall;[Score FOVres(iFOV).InfoTBL(:,2:end)]];
    TBLall_sub=[TBLall_sub;[Score_sub FOVres(iFOV).InfoTBL(:,2:end)]];

    
end

iWin=1;
writetable(TBLallNontarget,[SaveFunFOVDim 'NonTargetSLMScoreTrialDynWin' num2str(iWin) '.csv']);
writetable(TBLall,[SaveFunFOVDim 'AllSLMScoreTrialDynWin' num2str(iWin) '.csv']);


for iFOV=1:length(IndexFOVNeed)
    % if iFOV==6
    %    continue;
    % end

    tblFOV=SLMGroupTableTrial(SLMGroupTableTrial.Session==iFOV,:);
    TrialID=unique(tblFOV.TrialID);
    ScoreNontarget=[];
    Score=[];
    Score_sub=[];
    DataVec=[];
    clear DataVecTBL
    for iTrial=1:size(FOVres(iFOV).DataVec,2)
        DataVec(:,iTrial)= tblFOV.Response(tblFOV.TrialID==TrialID(iTrial),:);
        I1=find(tblFOV.TrialID==TrialID(iTrial));
    
        % I1=find(tblFOV.TrialID==TrialID(1));
        NonTarget=find(tblFOV.TargetCell(I1)==0);
        Target=find(tblFOV.TargetCell(I1)==1);
    
        % temp=tblFOV(tblFOV.TrialID==TrialID(iTrial),:);
        % DataVecInfo(iTrial,:)=[tblFOV.Group(I1(1)) tblFOV.Whisk(I1(1)) tblFOV.PowerZero(I1(1)) tblFOV.Time(I1(1))  tblFOV.Speed(I1(1))];
        % DataVec(tblFOV.TargetCell(I1)==1,iTrial)=0;
        DataVecTBL(iTrial,:)=tblFOV(I1(1),:);
        DataVecTBL(iTrial,:).SpeedR=mean(tblFOV.SpeedR(I1(NonTarget),:));
        DataVecTBL(iTrial,:).StimR=mean(tblFOV.StimR(I1(NonTarget),:));
        % DataVecTBL(iTrial,:).TargetSpeedR=mean(tblFOV.SpeedR(I1(Target),:));
        % DataVecTBL(iTrial,:).TargetStimR=mean(tblFOV.StimR(I1(Target),:));
    
        % 
        % %Keep Score positive
        Score(:,iTrial)=lsqnonneg(FOVres(iFOV).W(:,FOVres(iFOV).ComI), FOVres(iFOV).DataVec(:,iTrial));
        ScoreNontarget(:,iTrial)=lsqnonneg(FOVres(iFOV).WNontarget(:,FOVres(iFOV).ComI), FOVres(iFOV).DataVecNonTarget(:,iTrial));
        Score_sub(:,iTrial)=lsqnonneg(FOVres(iFOV).W_sub(:,FOVres(iFOV).ComI_sub), FOVres(iFOV).DataVecNonTarget(:,iTrial));
        % 
        % % %%Keep Score positive

    end
    Score=Score';
    ScoreNontarget=ScoreNontarget';
    Score_sub=Score_sub';

    % Score=zscore(Score);
    % ScoreNontarget=zscore(ScoreNontarget);
    % Score_sub=zscore(Score_sub);

    pinvW = pinv(FOVres(iFOV).W(:,FOVres(iFOV).ComI));
    % Score = (pinvW * FOVres(iFOV).DataVec)'; 
    % Score(:,1) = Score(:,1) - [ones(length(Score(:,1)),1) DataVecTBL.Speed]*FOVres(iFOV).SpeedReg.B;

% x=[min(data1(:)) max(data1(:))]';
% y=[[1 1]' x]*B;
% close all
% plot(lags,c_sub(:,12))
    subNonTargetData=FOVres(iFOV).DataVec;
    subNonTargetData(Target,:)=0;
    % ScoreNontarget = (pinvW * subNonTargetData)'; 
    % ScoreNontarget(:,1) =  ScoreNontarget(:,1) - [ones(length(Score(:,1)),1) DataVecTBL.Speed]*FOVres(iFOV).SpeedRegNontarget.B;

    pinvW_sub = pinv(FOVres(iFOV).W_sub(:,FOVres(iFOV).ComI_sub));
    % Score_sub = (pinvW_sub * FOVres(iFOV).DataVecNonTarget)'; 

    % Score_sub(:,1) = Score_sub(:,1) - [ones(length(Score(:,1)),1) DataVecTBL.Speed]*FOVres(iFOV).SpeedReg_sub.B;

    % ScoreNontarget=array2table(ScoreNontarget,"VariableNames",{'SpeedScore','WhiskerScore'});
    % Score=array2table(Score,"VariableNames",{'SpeedScore','WhiskerScore'});

    % ScoreNontarget=zscore(FOVres(iFOV).DataVecNonTarget'*FOVres(iFOV).WNontarget(:,FOVres(iFOV).ComI));
    % ScoreNontarget=array2table(ScoreNontarget',"VariableNames",{'SpeedScore','WhiskerScore'});
    % Score=array2table(Score',"VariableNames",{'SpeedScore','WhiskerScore'});
    % 
    ScoreNontarget=array2table(ScoreNontarget,"VariableNames",{'SpeedScore','WhiskerScore'});
    Score=array2table(Score,"VariableNames",{'SpeedScore','WhiskerScore'});
    Score_sub=array2table(Score_sub,"VariableNames",{'SpeedScore','WhiskerScore'});


    % ScoreNontarget=array2table(zscore(ScoreNontarget),"VariableNames",{'SpeedScore','WhiskerScore'});
    % Score=array2table(zscore(Score),"VariableNames",{'SpeedScore','WhiskerScore'});
    % Score_sub=array2table(zscore(Score_sub),"VariableNames",{'SpeedScore','WhiskerScore'});

    % ScoreNontarget=array2table(ScoreNontarget,"VariableNames",{'SpeedScore','WhiskerScore'});
    % Score=array2table(Score,"VariableNames",{'SpeedScore','WhiskerScore'});

    TBLallNontarget=[TBLallNontarget;[ScoreNontarget FOVres(iFOV).InfoTBL(:,2:end)]];
    TBLall=[TBLall;[Score FOVres(iFOV).InfoTBL(:,2:end)]];
    TBLall_sub=[TBLall_sub;[Score_sub FOVres(iFOV).InfoTBL(:,2:end)]];

    
end



   GroupPair.CorrName='fdr';
   GroupPair.Q=0.1;
   GroupPair.Pair=[];
   GroupPair.SignY=0.006;
   GroupPair.Plot=1;
   GroupPair.Std=0;      %%%%%%%%%using standard deviation as errorbar
   GroupPair.SamplePlot=1; %%%%%%%%%Plot Individual Sample Point
   GroupPair.SamplePairedPlot=1; %%%%%%%%%Dash line for paired comparison sample
   GroupPair.LimY=[-0.006 GroupPair.SignY*1.2];
   GroupPair.Marker={'o'};   % Datatype=varargin{1};  %%%%%%Datatype=0, non-paired data; 1 for paired data


GroupMethod='mean';
SpeedTh=1000
TBLallavg = groupsummary(TBLall(TBLall.Whisk==0&TBLall.PowerZero==0&abs(TBLall.Speed)<SpeedTh,:), {'Session', 'Group'}, GroupMethod);
TBL_subavg= groupsummary(TBLall_sub(TBLall_sub.Whisk==0&TBLall_sub.PowerZero==0&abs(TBLall_sub.Speed)<SpeedTh,:),{'Session', 'Group'}, GroupMethod);
% 
% TBLallavg = groupsummary(TBLall(TBLall.Whisk==0&TBLall.PowerZero==0&TBLall.Group~=3&abs(TBLall.Speed)<SpeedTh,:), {'Session', 'Group'}, GroupMethod);
% TBL_subavg= groupsummary(TBLall_sub(TBLall.Whisk==0&TBLall.PowerZero==0&TBLall.Group~=3&abs(TBLall.Speed)<SpeedTh,:),{'Session', 'Group'}, GroupMethod);
% 
% SpeedTh=1000;
% TBLallavg = groupsummary(TBLall(TBLall.Whisk==0&TBLall.PowerZero==0&TBLall.Group~=3&TBLall.Speed<SpeedTh,:), {'Session', 'Group'}, GroupMethod);
% TBL_subavg= groupsummary(TBLall_sub(TBLall.Whisk==0&TBLall.PowerZero==0&TBLall.Group~=3&TBLall.Speed<SpeedTh,:),{'Session', 'Group'}, GroupMethod);
% 
oldNames = TBLallavg.Properties.VariableNames;
% Find which names start with 'mean_'
meanVars = startsWith(oldNames, [GroupMethod '_']);
% Remove 'mean_' prefix
newNames = oldNames;
newNames(meanVars) = erase(oldNames(meanVars), [GroupMethod '_']);
TBLallavg.Properties.VariableNames=newNames;
TBL_subavg.Properties.VariableNames=newNames;



% for igroup=1:3
%     tempData{igroup,1}=TBLallavg.mean_SpeedScore(TBLallavg.Group==igroup);
% end
% a=ErrorBarPlotLU(1:3,tempData,[],ProcessPar.GroupColor,2,1,[],GroupPair,[1 1 1])
% 

ProcessPar.GroupColor

   Param.Color=ProcessPar.GroupColor;
   Param.Marker='o';
   Param.MarkerSize=10;
   Param.Rtype='pearson';
   Param.xLim=[min(TBLallavg.TargetSpeedR) max(TBLallavg.TargetSpeedR)];
   Param.yLim=[min(TBLallavg.SpeedScore) max(TBLallavg.SpeedScore)];
   Param.xLabel=[];
   Param.yLabel=[];

figure;
[OutPut,r,p,h]=LuPairRegressPlot_Group_Cov(TBLallavg.TargetSpeedR,TBLallavg.SpeedScore,TBLallavg.Speed,TBLallavg.Group,Param)
  
figure;
[OutPut,r,p,h]=LuPairRegressPlot_Group_Cov(TBL_subavg.TargetSpeedR,TBL_subavg.SpeedScore,TBL_subavg.Speed,TBL_subavg.Group,Param)

SpeedTh=1;

close all
Param2=Param;
Param2.Color=Param.Color(1:2,:);
figure;
[OutPut,r,p,h]=LuPairRegressPlot_Group_Cov(TBLallavg.TargetSpeedScore,TBLallavg.SpeedScore,TBLallavg.Speed,TBLallavg.Group,Param2)
  
figure;
[OutPut,r,p,h]=LuPairRegressPlot_Group_Cov(TBLall.TargetSpeedR,TBLall.SpeedScore,TBLall.Speed,TBLall.Group,Param)
figure;
[OutPut,r,p,h]=LuPairRegressPlot_Group_Cov(TBLall.TargetSpeedR,TBLall.SpeedScore,TBLall.Speed,TBLall.Group,Param)



SpeedTh=2;
TBLallTemp=TBLall(TBLall.Whisk==0&TBLall.PowerZero==0&abs(TBLall.Speed)<SpeedTh,:);
TBLsubTemp=TBLall_sub(TBLall_sub.Whisk==0&TBLall_sub.PowerZero==0&abs(TBLall_sub.Speed)<SpeedTh,:);
figure;
[OutPut,r,p,h]=LuPairRegressPlot_Group_Cov(TBLallTemp.TargetSpeedR,TBLallTemp.SpeedScore,TBLallTemp.Speed,TBLallTemp.Group,Param)
figure;
[OutPut,r,p,h]=LuPairRegressPlot_Group_Cov(TBLsubTemp.TargetSpeedR,TBLsubTemp.SpeedScore,TBLsubTemp.Speed,TBLsubTemp.Group,Param)


SpeedTh=2;
TBLallTemp=TBLall(TBLall.Whisk==0&TBLall.PowerZero==0&TBLall.Speed<SpeedTh,:);
TBLsubTemp=TBLall_sub(TBLall_sub.Whisk==0&TBLall_sub.PowerZero==0&TBLall_sub.Speed<SpeedTh,:);
figure;
[OutPut,r,p,h]=LuPairRegressPlot_Group_Cov(TBLallTemp.TargetSpeedR,TBLallTemp.SpeedScore,TBLallTemp.Speed,TBLallTemp.Group,Param)
figure;
[OutPut,r,p,h]=LuPairRegressPlot_Group_Cov(TBLsubTemp.TargetSpeedR,TBLsubTemp.SpeedScore,TBLsubTemp.Speed,TBLsubTemp.Group,Param)


TBLallTemp=TBLall(TBLall.Whisk==1&TBLall.PowerZero==0&TBLall.Speed<SpeedTh,:);
TBLsubTemp=TBLall_sub(TBLall_sub.Whisk==1&TBLall_sub.PowerZero==0&TBLall_sub.Speed<SpeedTh,:);

figure;
[OutPut,r,p,h]=LuPairRegressPlot_Group_Cov(TBLallTemp.TargetStimR,TBLallTemp.WhiskerScore,TBLallTemp.Speed,TBLallTemp.Group,Param)
figure;
[OutPut,r,p,h]=LuPairRegressPlot_Group_Cov(TBLsubTemp.TargetStimR,TBLsubTemp.WhiskerScore,TBLsubTemp.Speed,TBLsubTemp.Group,Param)

figure;
[OutPut,r,p,h]=LuPairRegressPlot_Group_Cov(TBL_subavg.TargetSpeedScore,TBL_subavg.SpeedScore,TBL_subavg.Speed,TBL_subavg.Group,Param)


GroupMethod='median';

TBLallavg = groupsummary(TBLall(TBLall.Whisk==1&TBLall.Group<=3&TBLall.PowerZero==0&TBLall.Speed<SpeedTh,:), {'Session', 'Group'}, GroupMethod);
TBL_subavg= groupsummary(TBLall_sub(TBLall.Whisk==1&TBLall.Group<=3&TBLall.PowerZero==0&TBLall.Speed<SpeedTh,:),{'Session', 'Group'}, GroupMethod);
oldNames = TBLallavg.Properties.VariableNames;
% Find which names start with 'mean_'
meanVars = startsWith(oldNames, [GroupMethod '_']);
% Remove 'mean_' prefix
newNames = oldNames;
newNames(meanVars) = erase(oldNames(meanVars), [GroupMethod '_']);
TBLallavg.Properties.VariableNames=newNames;
TBL_subavg.Properties.VariableNames=newNames;



Param.Color=ProcessPar.GroupColor;
Param.Marker='o';
Param.MarkerSize=10;
Param.Rtype='pearson';
Param.xLim=[min(TBLallavg.TargetStimR) max(TBLallavg.TargetStimR)];
Param.yLim=[min(TBLallavg.WhiskerScore) max(TBLallavg.WhiskerScore)];
figure;
[OutPut,r,p,h]=LuPairRegressPlot_Group_Cov(TBLallavg.TargetStimR,TBLallavg.WhiskerScore,TBLallavg.Speed,TBLallavg.Group,Param)
figure;
[OutPut,r,p,h]=LuPairRegressPlot_Group_Cov(TBL_subavg.TargetStimR,TBL_subavg.WhiskerScore,TBL_subavg.Speed,TBL_subavg.Group,Param)






SpeedTh=1;
% A=TBLallNontarget(TBLallNontarget.Whisk==0&TBLallNontarget.PowerZero==0&TBLallNontarget.Speed<SpeedTh,:);
% B=TBLall(TBLall.Whisk==0&TBLall.PowerZero==0&TBLall.Speed<SpeedTh,:);
figure;
clear dataAll dataNonT dataAll_sub statsSpeed
for iWhisk=1:length(ProcessPar.VolOut)

tempNonT=TBLallNontarget(TBLallNontarget.Whisk==ProcessPar.VolOut(iWhisk)&TBLallNontarget.PowerZero==0&TBLallNontarget.Speed<SpeedTh,:);
tempAll=TBLall(TBLall.Whisk==ProcessPar.VolOut(iWhisk)&TBLall.PowerZero==0&TBLall.Speed<SpeedTh,:);
tempAll_sub=TBLall_sub(TBLall_sub.Whisk==ProcessPar.VolOut(iWhisk)&TBLall_sub.PowerZero==0&TBLall_sub.Speed<SpeedTh,:);
% 
% tempNonT=TBLallNontarget(TBLallNontarget.Whisk==ProcessPar.VolOut(iWhisk)&TBLallNontarget.PowerZero==0&abs(TBLallNontarget.Speed)<SpeedTh,:);
% tempAll=TBLall(TBLall.Whisk==ProcessPar.VolOut(iWhisk)&TBLall.PowerZero==0&abs(TBLall.Speed)<SpeedTh,:);
% tempAll_sub=TBLall_sub(TBLall_sub.Whisk==ProcessPar.VolOut(iWhisk)&TBLall_sub.PowerZero==0&abs(TBLall_sub.Speed<SpeedTh),:);


for iGroup=1:3
    dataAll{iGroup,iWhisk}=tempAll.SpeedScore(tempAll.Group==iGroup);
    dataNonT{iGroup,iWhisk}=tempNonT.SpeedScore(tempNonT.Group==iGroup);
    dataAll_sub{iGroup,iWhisk}=tempAll_sub.SpeedScore(tempAll_sub.Group==iGroup);

end

subplotLU(2,3,iWhisk,1)
% statsSpeed(iWhisk,1)=ErrorBoxPlotLU(1:3,dataAll(:,iWhisk),ProcessPar.GroupColor,[],GroupPair)
statsSpeed(iWhisk,1)=ErrorBarPlotLU(1:3,dataAll(:,iWhisk),[],ProcessPar.GroupColor,2,0,[],GroupPair)

set(gca,'ylim',GroupPair.LimY,'ytick',GroupPair.LimY);
subplotLU(2,3,iWhisk,2)
% statsSpeed(iWhisk,2)=ErrorBoxPlotLU(1:3,dataNonT(:,iWhisk),ProcessPar.GroupColor,[],GroupPair)
statsSpeed(iWhisk,2)=ErrorBarPlotLU(1:3,dataNonT(:,iWhisk),[],ProcessPar.GroupColor,2,0,[],GroupPair)

set(gca,'ylim',GroupPair.LimY,'ytick',GroupPair.LimY);
subplotLU(2,3,iWhisk,3)
% statsSpeed(iWhisk,3)=ErrorBoxPlotLU(1:3,dataAll_sub(:,iWhisk),ProcessPar.GroupColor,[],GroupPair)
statsSpeed(iWhisk,3)=ErrorBarPlotLU(1:3,dataAll_sub(:,iWhisk),[],ProcessPar.GroupColor,2,0,[],GroupPair)

set(gca,'ylim',GroupPair.LimY,'ytick',GroupPair.LimY);



end
statsSpeed(1,1)
statsSpeed(1,3)

figure;
clear dataAll dataNonT dataAll_sub
for iWhisk=1:length(ProcessPar.VolOut)

tempNonT=TBLallNontarget(TBLallNontarget.Whisk==ProcessPar.VolOut(iWhisk)&TBLallNontarget.PowerZero==0&TBLallNontarget.Speed<SpeedTh,:);
tempAll=TBLall(TBLall.Whisk==ProcessPar.VolOut(iWhisk)&TBLall.PowerZero==0&TBLall.Speed<SpeedTh,:);
tempAll_sub=TBLall_sub(TBLall_sub.Whisk==ProcessPar.VolOut(iWhisk)&TBLall_sub.PowerZero==0&TBLall_sub.Speed<SpeedTh,:);
% 
tempNonT=TBLallNontarget(TBLallNontarget.Whisk==ProcessPar.VolOut(iWhisk)&TBLallNontarget.PowerZero==0&abs(TBLallNontarget.Speed)<SpeedTh,:);
tempAll=TBLall(TBLall.Whisk==ProcessPar.VolOut(iWhisk)&TBLall.PowerZero==0&abs(TBLall.Speed)<SpeedTh,:);
tempAll_sub=TBLall_sub(TBLall_sub.Whisk==ProcessPar.VolOut(iWhisk)&TBLall_sub.PowerZero==0&abs(TBLall_sub.Speed<SpeedTh),:);


for iGroup=1:3
    dataAll{iGroup,iWhisk}=tempAll.SpeedScore(tempAll.Group==iGroup);
    dataNonT{iGroup,iWhisk}=tempNonT.SpeedScore(tempNonT.Group==iGroup);
    dataAll_sub{iGroup,iWhisk}=tempAll_sub.SpeedScore(tempAll_sub.Group==iGroup);

end

subplotLU(2,3,iWhisk,1)
statsSpeed(iWhisk,1)=ErrorBoxPlotLU(1:3,dataAll(:,iWhisk),ProcessPar.GroupColor,[],GroupPair)
set(gca,'ylim',GroupPair.LimY,'ytick',GroupPair.LimY);
subplotLU(2,3,iWhisk,2)
statsSpeed(iWhisk,2)=ErrorBoxPlotLU(1:3,dataNonT(:,iWhisk),ProcessPar.GroupColor,[],GroupPair)
set(gca,'ylim',GroupPair.LimY,'ytick',GroupPair.LimY);
subplotLU(2,3,iWhisk,3)
statsSpeed(iWhisk,3)=ErrorBoxPlotLU(1:3,dataAll_sub(:,iWhisk),ProcessPar.GroupColor,[],GroupPair)
set(gca,'ylim',GroupPair.LimY,'ytick',GroupPair.LimY);



end
statsSpeed(1,3)


   GroupPair.SignY=0.003;
   GroupPair.Plot=1;
   GroupPair.Std=0;      %%%%%%%%%using standard deviation as errorbar
   GroupPair.SamplePlot=1; %%%%%%%%%Plot Individual Sample Point
   GroupPair.SamplePairedPlot=0; %%%%%%%%%Dash line for paired comparison sample
   GroupPair.LimY=[-0.001 GroupPair.SignY*1.2];

figure;
clear dataAll dataNonT dataAll_sub
for iWhisk=1:length(ProcessPar.VolOut)

tempNonT=TBLallNontarget(TBLallNontarget.Whisk==ProcessPar.VolOut(iWhisk)&TBLallNontarget.PowerZero==0&TBLallNontarget.Speed<SpeedTh,:);
tempAll=TBLall(TBLall.Whisk==ProcessPar.VolOut(iWhisk)&TBLall.PowerZero==0&TBLall.Speed<SpeedTh,:);
tempAll_sub=TBLall_sub(TBLall_sub.Whisk==ProcessPar.VolOut(iWhisk)&TBLall_sub.PowerZero==0&TBLall_sub.Speed<SpeedTh,:);




for iGroup=1:3
    dataAll{iGroup,iWhisk}=tempAll.WhiskerScore(tempAll.Group==iGroup);
    dataNonT{iGroup,iWhisk}=tempNonT.WhiskerScore(tempNonT.Group==iGroup);
    dataAll_sub{iGroup,iWhisk}=tempAll_sub.WhiskerScore(tempAll_sub.Group==iGroup);

end

subplotLU(2,3,iWhisk,1)
statsWhisker(iWhisk,1)=ErrorBoxPlotLU(1:3,dataAll(:,iWhisk),ProcessPar.GroupColor,[],GroupPair)
set(gca,'ylim',[-0.001 0.005],'ytick',[-0.001:0.001:0.004]);

subplotLU(2,3,iWhisk,2)
statsWhisker(iWhisk,2)=ErrorBoxPlotLU(1:3,dataNonT(:,iWhisk),ProcessPar.GroupColor,[],GroupPair)
set(gca,'ylim',[-0.001 0.005],'ytick',[-0.001:0.001:0.004]);

subplotLU(2,3,iWhisk,3)
statsWhisker(iWhisk,3)=ErrorBoxPlotLU(1:3,dataAll_sub(:,iWhisk),ProcessPar.GroupColor,[],GroupPair)
set(gca,'ylim',[-0.001 0.005],'ytick',[-0.001:0.001:0.004]);


end



SpeedTh=1;
A=TBLallNontarget(TBLallNontarget.Whisk==1&TBLallNontarget.PowerZero==0&abs(TBLallNontarget.Speed)<SpeedTh,:);
B=TBLall(TBLall.Whisk==1&TBLall.PowerZero==0&abs(TBLall.Speed)<SpeedTh,:);

for iGroup=1:3
    tempData{iGroup}=A.WhiskerScore(A.Group==iGroup);
    tempData2{iGroup}=B.WhiskerScore(B.Group==iGroup);
end

   GroupPair.CorrName='fdr';
   GroupPair.Q=0.1;
   GroupPair.Pair=[];
   GroupPair.SignY=2;
   GroupPair.Plot=1;
   GroupPair.Std=0;      %%%%%%%%%using standard deviation as errorbar
   GroupPair.SamplePlot=1; %%%%%%%%%Plot Individual Sample Point
   GroupPair.SamplePairedPlot=0; %%%%%%%%%Dash line for paired comparison sample
   GroupPair.LimY=[-1 GroupPair.SignY*1.2];
   GroupPair.Marker={'o'};   % Datatype=varargin{1};  %%%%%%Datatype=0, non-paired data; 1 for paired data
   % SavePath=varargin{2};  %%%%%%Save the stats to txt file
   % GroupPair=varargin{3};
   % RGroupID=varargin{4};
   

% figure;
% stats=ErrorBoxPlotLU(1:3,tempData,ProcessPar.GroupColor,[])

figure;
subplot(1,2,1)
stats=ErrorBoxPlotLU(1:3,tempData2,ProcessPar.GroupColor,[])
subplot(1,2,2)
stats=ErrorBoxPlotLU(1:3,tempData,ProcessPar.GroupColor,[])


figure;
ErrorBarPlotLU(1:3,dataAll_sub(:,1),[],ProcessPar.GroupColor,2,0,[])

figure;
subplot(1,2,1)
stats=ErrorBarPlotLU(1:3,tempData2,[],ProcessPar.GroupColor,2,0,[],GroupPair)
subplot(1,2,2)
stats=ErrorBarPlotLU(1:3,tempData,[],ProcessPar.GroupColor,2,0,[],GroupPair)

stats=ErrorViolinHalf(1:3,tempData,ProcessPar.GroupColor,GroupPair)

x=1:3
GroupPair.CorrName='fdr';
GroupPair.Test='Ttest';
GroupPair.Q=0.1;
GroupPair.Pair=[];
GroupPair.SignY=2;
GroupPair.Plot=1;
GroupPair.Std=0;      %%%%%%%%%using standard deviation as errorbar
GroupPair.SamplePlot=1; %%%%%%%%%Plot Individual Sample Point
GroupPair.SamplePairedPlot=0; %%%%%%%%%Dash line for paired comparison sample
GroupPair.LimY=[0 GroupPair.SignY*1.2];
GroupPair.Marker={'o'};
GroupPair.ViolinLR=zeros(1,length(x));
GroupPairSpeed=GroupPair;
GroupPairSpeed.SignY=3;
GroupPairSpeed.LimY=[-1 GroupPairSpeed.SignY*1.2];

figure;
stats=ErrorViolinHalf(1:3,tempData2,ProcessPar.GroupColor,0,[],GroupPair,[1 2 3])

figure;
stats=ErrorViolinHalf(1:3,tempData,ProcessPar.GroupColor,0,[],GroupPair,[1 2 3])
set(gca,'ylim',[-3 3])


iData=1;
pcaN=20;
ScoreLabel={'SpeedScore','WhiskScore'}
DataTable=[];
for iFOV=1:length(IndexFOVNeed)
tblFOV=SLMGroupTableTrial(SLMGroupTableTrial.Session==iFOV,:);
TrialID=unique(tblFOV.TrialID);
tbltemp=SLMGroupTableTrial(SLMGroupTableTrial.Session==iFOV&SLMGroupTableTrial.TrialID==TrialID(1),:);
% [W,H] = nnmf(NeuroTrace{iFOV}{iData}(tbltemp.TargetCell==0,SponInd),pcaN);

I1=find(tblFOV.TrialID==TrialID(1));
NonTarget=find(tblFOV.TargetCell(I1)==0);

StructToVars(FOVres(iFOV));
pcaIspeed=ComI(1);
pcaIstim=ComI(2);



tblFOV=SLMGroupTableTrial(SLMGroupTableTrial.Session==iFOV,:);
TrialID=unique(tblFOV.TrialID);
WSub=W(:,[pcaIspeed pcaIstim]);
WNontargetSub=WNontarget(:,[pcaIspeed pcaIstim]);




Score=DataVec'*WSub;
ScoreNontarget=DataVecNonTarget'*WNontargetSub;


InfoTable(:,2:end)




end






% close all
% for PCI=1:10
%     figure;
%     [~,rankI]=sort((W(:,PCI)),'descend');
%     imagesc(NeuroTrace{iFOV}{iData}(rankI,SponInd));
%     % subplot(1,2,1)
%     % plot(AmpNormalize(coeff(SponInd,PCI),[1 99]),'k');
%     % hold on;plot(AmpNormalize(nanmean(BehTrace(iFOV).Speed(SponInd,:),2)),'g')
%     % % plot(AmpNormalize(coeff(SponInd,1),[1 99]),'r');
%     % hold on;plot(AmpNormalize(nanmean(BehTrace(iFOV).Stim(SponInd,:),2)),'r')
% end




close all
% figure;
for PCI=6:11
    figure;
    % subplot(1,2,1)
    plot(AmpNormalize(H(SponInd,PCI),[1 99]),'k');
    hold on;plot(AmpNormalize(nanmean(BehTrace(iFOV).Speed(SponInd,:),2)),'g')
    % plot(AmpNormalize(coeff(SponInd,1),[1 99]),'r');
    hold on;plot(AmpNormalize(nanmean(BehTrace(iFOV).Stim(SponInd,:),2)),'r')

end


figure;
plot(lags,c(:,10))
figure;
plot(lags,c(:,2))






close all
for PCI=1:5
    figure;
    [~,rankI]=sort((score(:,PCI)),'descend');
    imagesc(NeuroTrace{iFOV}{iData}(rankI,SponInd));
    % subplot(1,2,1)
    % plot(AmpNormalize(coeff(SponInd,PCI),[1 99]),'k');
    % hold on;plot(AmpNormalize(nanmean(BehTrace(iFOV).Speed(SponInd,:),2)),'g')
    % % plot(AmpNormalize(coeff(SponInd,1),[1 99]),'r');
    % hold on;plot(AmpNormalize(nanmean(BehTrace(iFOV).Stim(SponInd,:),2)),'r')
end
[coeff,score,latent,tsquared,explained,mu] = pca(NeuroTrace{iFOV}{iData}(:,SponInd));
explained(1:10);


pcaN=10;
rSpeed=corr(coeff(ExcludeStimInd,1:pcaN),nanmean(BehTrace(iFOV).Speed(ExcludeStimInd,:),2),'rows','complete')

maxLag=10;
clear rStim c
for PCI = 1:pcaN
    [c(:,PCI), lags] = xcorr(coeff(SponInd,PCI), AmpNormalize(nanmean(BehTrace(iFOV).Stim(SponInd,:),2),[0 100]), maxLag, 'coeff');
    PostI = find(lags >= 0);
    PreI=find(lags<0);
    [~, i1] = max(abs(c(PostI,PCI))-mean(c(PreI,PCI)));
    rStim(PCI, 1) = c(PostI(i1),PCI);
end

[~,pcaIstim]=max(rStim)

[rankSpeedI,pcaIstim]=sort(rSpeed,'descend')

close all;
for i=1:12
figure;
plot(lags,c(:,i))
end
figure;
% plot(lags,c(:,1))
plot(lags,c(:,:))

close all
% figure;
for PCI=1:11
    figure;
    % subplot(1,2,1)
    plot(AmpNormalize(H(SponInd,PCI),[1 99]),'k');
    hold on;plot(AmpNormalize(nanmean(BehTrace(iFOV).Speed(SponInd,:),2)),'g')
    % plot(AmpNormalize(coeff(SponInd,1),[1 99]),'r');
    hold on;plot(AmpNormalize(nanmean(BehTrace(iFOV).Stim(SponInd,:),2)),'r')

end
























