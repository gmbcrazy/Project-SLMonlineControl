clear all
BatchSavePath='D:\Project1-LocalProcessing\Step1\';


load([BatchSavePath '03-Jul-2025FOV.mat'])
Suite2pDataKeywords='awakeRefSpon';

DataSavePath='\\nimhlabstore1.nimh.nih.gov\UFNC\FNC2\Zhang\Projects\Project-LocalProcessing\Step2\';
mkdir(DataSavePath);
DataSavePath=[DataSavePath Suite2pDataKeywords '\'];
mkdir(DataSavePath);
SaveFunCon=[DataSavePath 'FunCon\'];
mkdir(SaveFunCon)



PSTHparam.PreSLMCal = 10; 
PSTHparam.PostSLMCal = 15;
PSTHparam.pTh = 0.05; 
PSTHparam.TestMethod = 'ranksum';
PSTHparam.MPFrameJump = 2;
PSTHparam.TestStepFrame = 3;    %%post-slm frames for Test whether SLM works
PSTHparam.iData = 1;    %%post-slm frames for Test whether SLM works

%% Initial align behaviors with imaging, identify SLM target cells
for iFOV=1:length(FOVUpdate)
% iFOV=3;
    FOVtemp=FOVUpdate(iFOV);
    suite2pFOVPathLocalTemp=suite2pFOVPathLocal{iFOV};
    suite2pFOVPathLocalTemp=suite2pFOVPath{iFOV};
    
    OfflineCombineProcess_OneFOV(FOVtemp, Suite2pDataKeywords, suite2pFOVPathLocalTemp,PSTHparam);
end
GroupLabel={'L','S','N'};
nGroup=length(GroupLabel);
GroupColor=[255 51 153;91 20 212;121 247 111]/255;
% NodeColor=repmat([0.9 0.9 0.9],length(iscell),1);

%% Power test data
SaveFunCon=[DataSavePath 'FunCon\'];
mkdir(SaveFunCon)

for iFOV=12:length(FOVUpdate)
% iFOV=3;
    FOVtemp=FOVUpdate(iFOV);
    suite2pFOVPathLocalTemp=suite2pFOVPathLocal{iFOV};
    OfflinePowerTest_OneFOV(FOVtemp, Suite2pDataKeywords, suite2pFOVPathLocalTemp,PSTHparam,SaveFunCon)
    % OfflinePowerTest_OneFOVTempForSoohyun(FOVtemp, Suite2pDataKeywords, suite2pFOVPathLocalTemp,PSTHparam,SaveFunCon)

    % OfflinePowerTest_OneFOVTempForSoohyun
end




%% Group SLM data
SaveFunCon=[DataSavePath 'GroupSLM\'];
mkdir(SaveFunCon)
for iFOV=1:length(FOVUpdate)
    FOVtemp=FOVUpdate(iFOV);
    suite2pFOVPathLocalTemp=suite2pFOVPathLocal{iFOV};
    OfflineSLMAbsSpeedControl_OneFOV(FOVtemp, Suite2pDataKeywords, suite2pFOVPathLocalTemp,PSTHparam,SaveFunCon)
end

















GroupLabel={'L','S','N'};
GroupList=[1 2 3];

nGroup=length(GroupLabel);
GroupColor=[255 51 153;91 20 212;121 247 111]/255;
NodeColor=repmat([0.9 0.9 0.9],size(Output.NeuroPos3DMeta,1),1);
PowerZeroColor=[0.5 0.5 0.5];
PowerZero=[0 1];
PowerZeroLabel={'SLM','FakeSLM'};
VolOut=[0 1];
VolOutLabel={'NoWhisk','Whisk'}
AwakeState=[1];
AwakeStateLabel={'Awake'}
GroupMetaName=[GroupLabel {'FakeSLM'}];
GroupMetaColor=[GroupColor;PowerZeroColor];




TimBinFrame = -PSTHparam.PreSLMCal:PSTHparam.PostSLMCal-1;


RateParam.PlotType=3
RateParam.statisP=1;
RateParam.LegendShow=0;
RateParam.Legend=[];
RateParam.TimeRepeatAnova=1;
RateParam.GroupRepeatAnova = 0;
RateParam.Paired = 0;
RateParam.RepeatAnova= 1;
RateParam.BinName='Time';
RateParam.Bin=TimBinFrame+0.5;
RateParam.Ytick=[0 2 4];
RateParam.SigPlot='Anova';
RateParam.Q=0.1;

P.xLeft=0.1;        %%%%%%Left Margin
P.xRight=0.02;       %%%%%%Right Margin
P.yTop=0.06;         %%%%%%Top Margin
P.yBottom=0.1;      %%%%%%Bottom Margin
P.xInt=0.04;         %%%%%%Width-interval between subplots
P.yInt=0.04;         %%%%%%Height-interval between subplots


PostPreDiffSpeedTh=[0.5 1 2 10000];

SaveP1 = 'D:\Project1-LocalProcessing\Step3\awakeRefSpon\GroupSLM\'




PSTHparam.PreSLMCal=10;        %<----------------------------------------------------------------------------------Edit,Frame # before SLM to calculate baseline map
PSTHparam.PostSLMCal=3;        %<----------------------------------------------------------------------------------Edit,Frame # after SLM to calculate responsive map
% PSTHparam.YLim=[-50 600];       % %<-----------------------------------------For Suite2p based ROI signal only method
PSTHparam.pTh=0.05;             % %<-----------------------------------------For Suite2p based ROI signal only method
PSTHparam.TestMethod='ranksum'; % %<-----------------------------------------For Suite2p based ROI signal only method
PSTHparam.MPFrameJump=2; % %<-----------------------------------------For Suite2p based ROI signal only method
PSTHparam.iData=1; % %<-----------------------------------------DeltaF being used to calculate correlation or response

% for iFOV=1:length(FOVUpdate)
iFOV=1;
FOVtemp=FOVUpdate(iFOV);
suite2pFOVPathLocalTemp=suite2pFOVPathLocal{iFOV};
resultPaths = findAllFoldersKeyWords(suite2pFOVPathLocal{iFOV}, Suite2pDataKeywords);
OfflinePowerTest_OneFOV(FOVtemp, Suite2pDataKeywords, suite2pFOVPathLocalTemp,PSTHparam,SaveFunCon)



metaDataPath=Get_ExpDataFolder(resultPaths{1},'Step1Basic',{'Step1Meta.mat','BehAll.mat','.png'});
load([metaDataPath 'Step1Meta.mat'])
confSet=SLMPosInfo.confSetFinal

Zdepth = confSet.scan_Z + confSet.ETL;
suite2pPath = findAllFoldersKeyWords(suite2pFOVPathLocalTemp, Suite2pDataKeywords);

confSet.save_path0=suite2pPath{1};
[NeuronPos3D,NeuronPos3DRaw,CaData,CaDataPlane,Neuronstat]=Extract_Suite2p(confSet);
% load('C:\Users\zhangl33\Projects\Project-SLMonlineControl\CodeAndPackages\subfun\Color\colorMapPN3.mat','colorMapPN1')
[~, ~, CellMedCenter, cellBoundary, ~] = Suite2pCellIDMapFromStat(CaData.statCell, [confSet.SLM_Pixels_Y confSet.SLM_Pixels_X]);
Cell3DPos=[CellMedCenter Zdepth(CaData.CellPlaneID)'];
CellDist=squareform(pdist(Cell3DPos));

iData=1;
PSTHparam.TestStepFrame=3;
iscell=find(CaData.iscell(CaData.iscell(:,1)>0));
PVpower=xmlPower2PVpower(SLMTestInfo.confSet.UncagingLaserPower);
PVpower=intersect(round(PVpower),SLMInfoTable.UncagingLaserPower);

% [AlignedtempNData,AlignedInfoTable,statCellRes,TargetResponse,TargetCellResP,TargetCellResR,CellSampleN]=Aligned_FromSuite2p(NData{iData},TargetCellList,iscell,Suite2pTable,PVpower,PSTHparam);
[AlignedtempNData,AlignedInfoTable,CellResponse,statCellRes,TargetResponse,TargetCellResP,TargetCellResR,CellSampleN]=Aligned_FromSuite2p(NData{iData},TargetCellList,Suite2pTable,PVpower,PSTHparam);

% end
pAll=[];
for iCell=1:length(TargetCellList)
    if SuccTarget(iCell)
       temp=statCellRes(iCell,PowerTargetI(iCell)).p(:);
       % temp(TargetCellList(iCell))=[];
       pAll=[pAll;temp];
    end
end
FDR=0.1;
[h, crit_pAll, adj_p]=fdr_bh(pAll,FDR,'pdep');crit_pAll
ParamNet.crit_pAll=crit_pAll;
ParamNet.GroupColor=GroupColor;
ParamNet.TargetCellList=TargetCellList;
ParamNet.TargetCellListFunGroup=TargetCellListFunGroup;
ParamNet.iscell=iscell;
ParamNet.SuccTarget=SuccTarget;
ParamNet.ScoreLim=[-0.25 0.25];
ParamNet.ResponseLim=[-0.3 0.3];
ParamNet.ScoreMap=slanCM('wildfire',64);
ParamNet.ScoreLabel='Speed Corr';
ParamNet.ResponseMap=slanCM('seismic',64);
ParamNet.NodeColor=NodeColor;
ParamNet.PowerTargetI=PowerTargetI;
ParamNet.statCellRes=statCellRes;
ParamNet.crit_pAll=crit_pAll;
ParamNet.SuccAmp=SuccAmp;
ParamNet.xMat=[0.01 0.2 0.2 0.2 0.2];
ParamNet.yMat=[0.7 0.7 0.7 0.7 0.7];


PowerTestAdj=zeros(length(iscell))+NaN;
NodeTarget=zeros(length(iscell),1);
DistTarget=zeros(length(iscell));


% crit_pAll=0.001;
iCount=0;
for iCell=1:length(TargetCellList)
    if SuccTarget(iCell)
       TargetC=TargetCellList(iCell);
       temp1=statCellRes(iCell,PowerTargetI(iCell)).delta;
       temp2=statCellRes(iCell,PowerTargetI(iCell)).p>crit_pAll;
       % temp2=[];
       temp1(temp2)=NaN;
       PowerTestAdj(TargetCellList(iCell),:)=temp1;
       NodeTarget(TargetCellList(iCell))=SuccAmp(iCell);
       DistTarget(TargetCellList(iCell),:)=CellDist(TargetCellList(iCell),:);
       iCount=iCount+1;
    end
end
PowerTestAdj=PowerTestAdj-diag(diag(PowerTestAdj));
PowerTestNode=NodeTarget;



iCount=0;
CellN=size(NeuronPos3D,1);
TargetNameDep=[1:size(AlignedtempNData,1)];
TargetNameDep=[];
StimFunGDep=[];

DistTargetDep=[];
DistPlaneTargetDep=[];
DataResponse=[];
ResonseSpeedDep=[];
ResonseStimDep=[];
MixedData=[];
for iCell=1:length(TargetCellList)
    if SuccTarget(iCell)
       TargetC=TargetCellList(iCell);
       temp1=statCellRes(iCell,PowerTargetI(iCell)).delta;
       temp2=statCellRes(iCell,PowerTargetI(iCell)).p>crit_pAll;
       % temp2=[];
       temp1(temp2)=NaN;

       DataResponse=[DataResponse;temp1(:)];

       PowerTestAdj(TargetCellList(iCell),:)=temp1;
       NodeTarget(TargetCellList(iCell))=SuccAmp(iCell);
       DistTarget(TargetCellList(iCell),:)=CellDist(TargetCellList(iCell),:);
       iCount=iCount+1;

       TargetNameDep=[TargetNameDep;[1:CellN]'];
       StimFunGDep=[StimFunGDep;zeros(CellN,1)+TargetCellListFunGroup(iCell)];
       ResonseSpeedDep=[ResonseSpeedDep;squeeze(CorrResults.CorrResults.rSpeed(:,1,iData))];
       ResonseStimDep=[ResonseStimDep;squeeze(CorrResults.CorrResults.rStim(:,1,iData))];


       DistTargetDep = [DistTargetDep; CellDist(TargetCellList(iCell),:)'];
       
       tempDist=CellDist(TargetCellList(iCell),:);
       NonSamePlaneI=find(abs(NeuronPos3D(:,3)-NeuronPos3D(TargetCellList(iCell),3)))>0.1;
       tempDist(NonSamePlaneI)=NaN;
       DistPlaneTargetDep=[DistPlaneTargetDep;tempDist(:)];
    end
end


MixedData =[DataResponse DistTargetDep DistPlaneTargetDep StimFunGDep ResonseSpeedDep ResonseStimDep TargetNameDep];
MixedData(:,end+1) = iFOV;

tbl = table(DataResponse, DistTargetDep, DistPlaneTargetDep, StimFunGDep, ...
    ResonseSpeedDep, ResonseStimDep, TargetNameDep, MixedData(:,end), ...
    'VariableNames', {'Response', 'Dist','DistPlane','StimGroup','SpeedR','StimR','Cell','Session'});

% save([SaveFunCon 'tempMixedLinearM.mat'], 'tbl');
% 
% 
% Invalid=isnan(tbl.Response)|tbl.Dist<=0.1;

temptbl=tbl;
temptbl(Invalid,:)=[];


[~,Session,~]=fileparts(suite2pFOVPathLocal{iFOV}(1:end-1));

SaveP1=[SaveFunCon Session '\'];
mkdir(SaveP1);

% csvwrite([SaveP1 'tempMixedLinearM.csv'], 'tbl');
writetable(tbl,[SaveFunCon 'PowerTestResponse.csv']);
% save([SaveFunCon 'tempMixedLinearM.mat'], 'tbl');


% % RProgramfolder='C:\Users\zhangl33\OneDrive - National Institutes of Health\Documents\R\R-3.5.0\bin';
% % 
% % Rfolder = 'C:\Users\lzhang481\ToolboxAndScript\GenMatCode\Statistics\ANOVAandMixEffect\R\';
% % save([Rfolder 'R_MixAnova2F1W1B.mat'],'Data','Dependent','-V6');
% % 
% % 
% % RunRcode([Rfolder 'MixAnova2F1W1B.R'],'C:\Program Files\R\R-3.5.0\bin')
% % ResultFile='C:\Users\lzhang481\ToolboxAndScript\GenMatCode\Statistics\ANOVAandMixEffect\R\MixAnova2F1W1B.txt';
% % movefile(ResultFile,[SavePath '_MixAnova2F1W1B.txt']);


DistScoreTemp=sum(DistTarget)/sum(SuccTarget)*SLMPosInfo.yaml.umPerlPixelX;

ParamNetDist=ParamNet;
ParamNetDist.AllCellRes=1; 
ParamNetDist.ScoreLim=[0;450];
ParamNetDist.ScoreMap=slanCM('amethyst',64);



close all

figure;
ParamNetDist.AllCellRes=1; 
ParamNetDist.ScoreLabel='Dist. from Target (µm)';
[GgraphOut,Res,r,p]=ResPosNetwork_ScoreNodeOrder(SLMPosInfo,PowerTestAdj,PowerTestNode,ParamNetDist,DistScoreTemp(:))
GgraphOut.p.ArrowSize=4;
% GgraphOut.p.MarkerSize=2;
GgraphOut.p.LineWidth=0.5;
papersizePX=[0 0 33 9];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
LuFontStandard
saveas(gcf,[SaveP1 'SingleSLMVsSpeed'],'tif');
% saveas(gcf,[ResultFolder 'SingleSLMVsSpeed.eps'],'epsc');
% saveas(gcf,[ResultFolder 'SingleSLMVsSpeed'],'fig');



ParamNetDist.AllCellRes=0; 
for iFun=1:length(GroupLabel)

    tempAdj=PowerTestAdj;
    tempAdj(TargetCellList(TargetCellListFunGroup~=iFun),:)=NaN;

    tempNode=PowerTestNode;
    tempNode(TargetCellList(TargetCellListFunGroup~=iFun))=0;
    [GgraphOut,Res,r,p]=ResPosNetwork_ScoreNodeOrder(SLMPosInfo,tempAdj,tempNode,ParamNetDist,DistScoreTemp(:))
GgraphOut.p.ArrowSize=4;
% GgraphOut.p.MarkerSize=2;
GgraphOut.p.LineWidth=0.5;

papersizePX=[0 0 33 9];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
LuFontStandard
% print(gcf, [SaveP1 GroupLabel{iFun} 'SingleSLMVsDist.svg'], '-dsvg', '-painters');
% print(gcf, [SaveP1 GroupLabel{iFun} 'SingleSLMVsSpeed.eps'], '-depsc', '-painters');
% print(gcf, [SaveP1 GroupLabel{iFun} 'SingleSLMVsDist.tif'], '-dtiffn', '-painters');
% print(gcf, [SaveP1 'SingleSLMVsSpeed.png'], '-dpng', '-painters');
% 
% 
% close all

end




rSpeed=CorrResults.rSpeed;
rStim=CorrResults.rStim;



ParamNet.AllCellRes=1; 
close all
[GgraphOut,Res,r,p]=ResPosNetwork_ScoreNodeOrder(SLMPosInfo,PowerTestAdj,PowerTestNode,ParamNet,rSpeed(:,1,iData))
GgraphOut.p.ArrowSize=4;
% GgraphOut.p.MarkerSize=2;
GgraphOut.p.LineWidth=0.5;
papersizePX=[0 0 33 9];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
% saveas(gcf,[SaveP1 'SingleSLMVsSpeed'],'tif');
% saveas(gcf,[SaveP1 'SingleSLMVsSpeed.eps'],'epsc');
% saveas(gcf,[SaveP1 'SingleSLMVsSpeed'],'fig');
LuFontStandard
print(gcf, [SaveP1 'SingleSLMVsSpeed.svg'], '-dsvg', '-painters');
print(gcf, [SaveP1 'SingleSLMVsSpeed.eps'], '-depsc', '-painters');
print(gcf, [SaveP1 'SingleSLMVsSpeed.tif'], '-dtiffn', '-painters');
% print(gcf, [SaveP1 'SingleSLMVsSpeed.png'], '-dpng', '-painters');

ParamNet.AllCellRes=0; 
for iFun=1:length(GroupLabel)

    tempAdj=PowerTestAdj;
    tempAdj(TargetCellList(TargetCellListFunGroup~=iFun),:)=NaN;

    tempNode=PowerTestNode;
    tempNode(TargetCellList(TargetCellListFunGroup~=iFun))=0;
    [GgraphOut,Res,r,p]=ResPosNetwork_ScoreNodeOrder(SLMPosInfo,tempAdj,tempNode,ParamNet,rSpeed(:,1,iData))
GgraphOut.p.ArrowSize=4;
% GgraphOut.p.MarkerSize=2;
GgraphOut.p.LineWidth=0.5;

papersizePX=[0 0 33 9];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
LuFontStandard
% print(gcf, [SaveP1 GroupLabel{iFun} 'SingleSLMVsSpeed.svg'], '-dsvg', '-painters');
% print(gcf, [SaveP1 GroupLabel{iFun} 'SingleSLMVsSpeed.eps'], '-depsc', '-painters');
% print(gcf, [SaveP1 GroupLabel{iFun} 'SingleSLMVsSpeed.tif'], '-dtiffn', '-painters');
% print(gcf, [SaveP1 'SingleSLMVsSpeed.png'], '-dpng', '-painters');
% 
% 
% close all

end





ParamNet.AllCellRes=0; 
P.xLeft=0.1;
P.yBottom=0.12;
P.yTop=0.04;
P.xInt=0.04;
P.yInt=0.1;

ParamNet.AllCellRes=1; 

ParamNet.ScoreLim=[-0.1 0.1];
ParamNet.ScoreLabel='Stim Corr.';
[GgraphOut,Res,r,p]=ResPosNetwork_ScoreNodeOrder(SLMPosInfo,PowerTestAdj,PowerTestNode,ParamNet,rStim(:,1,iData))
GgraphOut.p.ArrowSize=8;
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
print(gcf, [SaveP1 'SingleSLMVsStim.svg'], '-dsvg', '-painters');
print(gcf, [SaveP1 'SingleSLMVsStim.eps'], '-depsc', '-painters');
print(gcf, [SaveP1 'SingleSLMVsStim.tif'], '-dtiff', '-painters');
papersizePX=[0 0 22 16];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
% LuFontStandard
% print(gcf, [SaveP1 'SepFunSingleSLMVsSpeedRegress.svg'], '-dsvg', '-painters');
% print(gcf, [SaveP1 'SepFunSingleSLMVsSpeedRegress.eps'], '-depsc', '-painters');
% print(gcf, [SaveP1 'SepFunSingleSLMVsSpeedRegress.tif'], '-dtiffn', '-painters');
% % print(gcf, [SaveP1 'SingleSLMVsSpeed.png'], '-dpng', '-painters');



close all
ParamNet.AllCellRes=0; 
for iFun=1:length(GroupLabel)

    tempAdj=PowerTestAdj;
    tempAdj(TargetCellList(TargetCellListFunGroup~=iFun),:)=NaN;

    tempNode=PowerTestNode;
    tempNode(TargetCellList(TargetCellListFunGroup~=iFun))=0;
    [GgraphOut,Res,r,p]=ResPosNetwork_ScoreNodeOrder(SLMPosInfo,tempAdj,tempNode,ParamNet,rStim(:,1,iData))
GgraphOut.p.ArrowSize=4;
% GgraphOut.p.MarkerSize=2;
GgraphOut.p.LineWidth=0.5;

papersizePX=[0 0 33 9];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
LuFontStandard
% print(gcf, [SaveP1 GroupLabel{iFun} 'SingleSLMVsStim.svg'], '-dsvg', '-painters');
% print(gcf, [SaveP1 GroupLabel{iFun} 'SingleSLMVsSpeed.eps'], '-depsc', '-painters');
% print(gcf, [SaveP1 GroupLabel{iFun} 'SingleSLMVStim.tif'], '-dtiffn', '-painters');
% print(gcf, [SaveP1 'SingleSLMVsSpeed.png'], '-dpng', '-painters');
% 
% 
% close all

end







ParamNet.AllCellRes=0; 

figure;
for iFun=1:length(GroupLabel)

    tempAdj=PowerTestAdj;
    tempAdj(TargetCellList(TargetCellListFunGroup~=iFun),:)=0;

    PosRate=nansum(tempAdj.*(tempAdj>0))';
    NegRate=nansum(tempAdj.*(tempAdj<0))';

   ResParam.Color=[0.1 0.1 0.1];
   ResParam.Marker='o';
   ResParam.MarkerSize=6;
   ResParam.Rtype='spearman';
   ResParam.xLim=ParamNet.ScoreLim;
   ResParam.yLim=[0 ParamNet.ResponseLim(2)];
   ResParam.xLabel=ParamNet.ScoreLabel;

   if iFun==1
   ResParam.yLabel='Excitation';
   else
   ResParam.yLabel='';
   end


subplotLU(2,length(GroupLabel),1,iFun,P);
if ParamNet.AllCellRes
[Res,r,p]=LuPairRegressPlot(squeeze(rStim(:,1,2)),PosRate(:),ResParam);
else
[Res,r,p]=LuPairRegressPlot(squeeze(rStim(PosRate>0,1,2)),PosRate(PosRate>0),ResParam);
end
   if iFun~=1
   set(gca,'yticklabel',[])
   end
   set(gca,'xticklabel',[])
   xlabel('')


   if iFun==1
   ResParam.yLabel='Inhibition';
   else
   ResParam.yLabel='';
   end

% ResParam.yLabel='Inhibition';
ResParam.yLim=[ParamNet.ResponseLim(1) 0];

subplotLU(2,length(GroupLabel),2,iFun,P);
if ParamNet.AllCellRes
[Res,r,p]=LuPairRegressPlot(squeeze(rStim(:,1,iData)),NegRate(:),ResParam);
else
[Res,r,p]=LuPairRegressPlot(squeeze(rStim(NegRate<0,1,2)),NegRate(NegRate<0),ResParam);
end
   if iFun~=1
   set(gca,'yticklabel',[])
   end
end
papersizePX=[0 0 22 16];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
LuFontStandard
% print(gcf, [SaveP1 'SepFunSingleSLMVsStimRegress.svg'], '-dsvg', '-painters');
% print(gcf, [SaveP1 'SepFunSingleSLMVsStimRegress.eps'], '-depsc', '-painters');
% print(gcf, [SaveP1 'SepFunSingleSLMVsStimRegress.tif'], '-dtiffn', '-painters');
% % print(gcf, [SaveP1 'SingleSLMVsSpeed.png'], '-dpng', '-painters');







% % % % % 
% % % % % 
% % % % % Suite2pDataKeywords, suite2pFOVPathLocal
% % % % % 
% % % % % for iFOV=1:length(FOVUpdate)
% % % % % FOVtemp=FOVUpdate(iFOV);
% % % % % 
% % % % % DataFolder=FOVtemp.DataFolder;
% % % % % ProcessFolder = [fileparts(FOVtemp.DataFolder(1:end-1)) '\']
% % % % % Temp1 = [fileparts(ProcessFolder(1:end-1)) '\']
% % % % % WorkFolder = [fileparts(Temp1(1:end-1)) '\']
% % % % % 
% % % % % 
% % % % % 
% % % % % SLMPosInfo=load([ProcessFolder 'SLMFunGroup.mat']);
% % % % % SLMTestInfo=load([ProcessFolder 'SLMIncludedIndFromIscell.mat']);
% % % % % confSet=SLMPosInfo.confSetFinal;
% % % % % Zdepth=confSet.scan_Z+confSet.ETL;
% % % % % Pos3DFun=SLMPosInfo.FinalPos3D;
% % % % % % FunScore=SLMPosInfo.FinalFunScore;
% % % % % Group=SLMPosInfo.Group;
% % % % % for iGroup=1:length(Group)
% % % % %     Pos3DGroup{iGroup}=Pos3DFun(Group(iGroup).Indices,:);
% % % % % end
% % % % % confSet.save_path0=DataFolder;
% % % % % Zdepth=confSet.scan_Z+confSet.ETL;
% % % % % nPlane=length(Zdepth);
% % % % % 
% % % % % Suite2pTable=FOVtemp.Suite2pTable;
% % % % % Suite2pTable=Suite2pTable(Suite2pTable.AwakeState<2,:);        %%%%%%%%%%Exclude anesia data
% % % % % 
% % % % % 
% % % % % 
% % % % % [FileIDList,i1]=unique(Suite2pTable.FileID);
% % % % % for i=1:length(FileIDList)
% % % % %     [fSpeed{i},fStim{i},timeStampCa{i},FrameTS{i},fVideo{i},VideoStartFrameTime{i}]=PV_VolExtract_MultiCyc(confSet,FileIDList(i));
% % % % % end
% % % % % 
% % % % % RemoveFrame=2;
% % % % % 
% % % % % SessFileTable=Suite2pTable(i1,:);
% % % % % tempTiffmark=[];
% % % % % suite2pInd=[];
% % % % % 
% % % % % CumTiffNum=cumsum(SessFileTable.Suite2pTiffNum);
% % % % % CumFrameNum=[0;CumTiffNum(1:end-1)]/nPlane;
% % % % % 
% % % % % for i=1:length(FileIDList)
% % % % %     i2=find(Suite2pTable.FileID==FileIDList(i));
% % % % %     for j=1:length(i2)
% % % % %         AddTemp=Suite2pTable.markCycle(i2(j))-(j-1)*RemoveFrame;
% % % % %         tempTiffmark=[tempTiffmark;AddTemp];
% % % % %         suite2pInd=[suite2pInd;AddTemp+CumFrameNum(i)];
% % % % %     end
% % % % % end
% % % % % clear i2
% % % % % 
% % % % % Suite2pTable.PostSLMTiffCycle=tempTiffmark;
% % % % % Suite2pTable.suite2pInd=suite2pInd;
% % % % % 
% % % % % 
% % % % % 
% % % % % resultPaths = findAllFoldersKeyWords(suite2pFOVPathLocal{iFOV}, Suite2pDataKeywords)
% % % % % if length(resultPaths)>1
% % % % %    disp('more than one suite2p data folder!')
% % % % % elseif length(resultPaths)==1
% % % % %    suite2pDataFolder=resultPaths{1};
% % % % % else
% % % % % 
% % % % % end
% % % % % confSet.save_path0=suite2pDataFolder
% % % % % 
% % % % % ResultFolder=[suite2pDataFolder 'result\'];
% % % % % mkdir(ResultFolder)
% % % % % ResultFolder=[ResultFolder 'Step1Basic\'];
% % % % % mkdir(ResultFolder)
% % % % % 
% % % % % 
% % % % % [NeuronPos3D,NeuronPos3DRaw,CaData,CaDataPlane,Neuronstat]=Extract_Suite2p(confSet);
% % % % % 
% % % % % TempconfSet=confSet;
% % % % % TempconfSet.save_path0=WorkFolder;
% % % % % [~,~,OnlineCaData,~,~]=Extract_Suite2p(TempconfSet);
% % % % % [~, ~, ~, cellBoundary, ~] = Suite2pCellIDMapFromStat(CaData.statCell, [512 512]);
% % % % % 
% % % % % SpeedAll=[];
% % % % % StimAll=[];
% % % % % for i=1:length(FileIDList)
% % % % %     i2=find(Suite2pTable.FileID==FileIDList(i));
% % % % % 
% % % % %     frameN=SessFileTable.Suite2pTiffNum(i)/nPlane;
% % % % %     if ~isempty(fSpeed{i})
% % % % %         InValid=[];
% % % % %         for j=1:RemoveFrame  
% % % % %             InValid=[InValid;Suite2pTable.markCycle(i2)+j-1];
% % % % %         end
% % % % %         InValid(isnan(InValid))=[];
% % % % %         tempSpeed=fSpeed{i};
% % % % %         tempStim=fStim{i};
% % % % %         tempSpeed(InValid,:)=[];
% % % % %         tempStim(InValid+1,:)=[];
% % % % %     else
% % % % %         tempSpeed=zeros(frameN,nPlane)+NaN;
% % % % %         tempStim=tempSpeed;
% % % % % 
% % % % %     end
% % % % %     SpeedAll=[SpeedAll;tempSpeed];
% % % % %     StimAll=[StimAll;tempStim];
% % % % % 
% % % % %     clear tempSpeed tempStim
% % % % % end
% % % % % 
% % % % % 
% % % % % 
% % % % % if size(SpeedAll,1) == size(CaData.F,2)
% % % % %    disp(['Speed length match neural data length'])
% % % % % else
% % % % %    disp(['Warning !! Speed length does Not match neural data length'])
% % % % % 
% % % % % end
% % % % % 
% % % % % SLMPos3D=SLMTestInfo.Pos3Dneed;
% % % % % SLMGroup=SLMTestInfo.FunScore(:,1);
% % % % % 
% % % % % 
% % % % % % DistTh=10;
% % % % % % [SLMtarget,SLMtargetCellDist]=SLMtargetMatchCell(SLMPos3D,NeuronPos3D,DistTh)
% % % % % % PointList1=find(SLMtarget>0);
% % % % % % 
% % % % % % 
% % % % % % RealTargetI=~isnan(SLMPosInfo.FinalFunScore(:,2))
% % % % % % SLMPos3D=SLMPosInfo.FinalPos3D(RealTargetI,:);
% % % % % % SLMGroup=SLMPosInfo.FinalFunScore(RealTargetI,1);
% % % % % % [SLMtarget,SLMtargetCellDist]=SLMtargetMatchCell(SLMPos3D,NeuronPos3D,DistTh)
% % % % % 
% % % % % 
% % % % % 
% % % % % Img=permute(CaData.PlaneMeanImg,[2,1,3]);
% % % % % ImgOnline=permute(OnlineCaData.PlaneMeanImgE,[2,1,3]);
% % % % % 
% % % % % % figure;
% % % % % % for iplane = 1:nPlane
% % % % % %     subplot(1,3,iplane)
% % % % % %     imshowpair(CaData.PlaneMeanImg(:,:,iplane),OnlineCaData.PlaneMeanImg(:,:,iplane))
% % % % % %     tform = imregcorr(CaData.PlaneMeanImgE(:,:,iplane),OnlineCaData.PlaneMeanImgE(:,:,iplane))
% % % % % % end
% % % % % 
% % % % % 
% % % % % % validtarget=SLMtarget>0;
% % % % % % figure;
% % % % % % H=MultiPlanes2DShow(Img, cellBoundary(SLMtarget(validtarget)), SLMPos3D(validtarget,:), round(SLMtargetCellDist(validtarget)), Zdepth, [1 0 0], [0 800])
% % % % % % figure;
% % % % % % H=MultiPlanes2DShow(ImgOnline, cellBoundary(SLMtarget), SLMPos3D, round(SLMtargetCellDist), Zdepth, [1 0 0], [0 800])
% % % % % 
% % % % % DistTh=6;
% % % % % [SLMtarget,SLMtargetCellDist]=SLMtargetMatchCell(SLMPos3D,NeuronPos3D,DistTh)
% % % % % PointList1=find(SLMtarget>0);
% % % % % 
% % % % % SLMtargetTable=zeros(size(Suite2pTable,1),1);
% % % % % SLMtargetTableGroup=zeros(size(Suite2pTable,1),1);
% % % % % 
% % % % % for iTarget=1:length(SLMtarget)
% % % % %     SLMtargetTable(find(Suite2pTable.Point==iTarget))=SLMtarget(iTarget);
% % % % %     SLMtargetTableGroup(find(Suite2pTable.Point==iTarget))=SLMGroup(iTarget);
% % % % % end
% % % % % Suite2pTable.PointTargetCell=SLMtargetTable;
% % % % % Suite2pTable.PointTargetCellGroup=SLMtargetTableGroup;
% % % % % 
% % % % % 
% % % % % [SLMFinalInSLMtest,~]=SLMtargetMatchCell(SLMPosInfo.FinalPos3D,SLMPos3D,0.1);
% % % % % % [SLMPosInfo.FinalPos3D(1:2,:) SLMPos3D(SLMFinalInSLMtest(1:2),:)]
% % % % % 
% % % % % 
% % % % % [TargetCellList,ia]=unique(SLMtargetTable(SLMtargetTable>0));
% % % % % temp=SLMtargetTableGroup(SLMtargetTable>0);
% % % % % TargetCellListFunGroup=temp(ia);
% % % % % 
% % % % % GroupTargetCell={};
% % % % % for iFun=1:length(SLMPosInfo.Group)
% % % % %     [temp,~]=SLMtargetMatchCell(SLMPosInfo.FinalPos3D(SLMPosInfo.Group(iFun).Indices,:),NeuronPos3D,DistTh);
% % % % %     GroupTargetCell{iFun}=temp(temp>0);
% % % % % end
% % % % % 
% % % % % 
% % % % % deltaFoF=double(F2deltaFoF(CaData.F,CaData.Fneu,CaData.ops.fs));
% % % % % iscell=find(CaData.iscell(:,1)>0.99);
% % % % % deltaFoF=deltaFoF(:,iscell);
% % % % % spks=double(CaData.spks(iscell,:));
% % % % % 
% % % % % 
% % % % % PSTHparam.PreSLMCal=10;        %<----------------------------------------------------------------------------------Edit,Frame # before SLM to calculate baseline map
% % % % % PSTHparam.PostSLMCal=3;        %<----------------------------------------------------------------------------------Edit,Frame # after SLM to calculate responsive map
% % % % % % PSTHparam.YLim=[-50 600];       % %<-----------------------------------------For Suite2p based ROI signal only method
% % % % % PSTHparam.pTh=0.05;             % %<-----------------------------------------For Suite2p based ROI signal only method
% % % % % PSTHparam.TestMethod='ranksum'; % %<-----------------------------------------For Suite2p based ROI signal only method
% % % % % PSTHparam.MPFrameJump=2; % %<-----------------------------------------For Suite2p based ROI signal only method
% % % % % 
% % % % % 
% % % % % TimBinFrame=-PSTHparam.PreSLMCal:PSTHparam.PostSLMCal-1;        %<----------------------------------------------------------------------------------Edit,Frame # before SLM to calculate baseline map
% % % % % 
% % % % % 
% % % % % NData={deltaFoF' spks};
% % % % % Nlabel={'DeltaF','Spks'};
% % % % % 
% % % % % 
% % % % % % CorrResults.rSpeedInitial=corr(deltaFoF,nanmean(SpeedAll,2),'rows','complete');
% % % % % % CorrResults.rSpeedInitial=corr(deltaFoF,nanmean(SpeedAll,2),'rows','complete');
% % % % % 
% % % % % Suite2pTable.Suite2pTiffNum(isnan(Suite2pTable.markCycle))/nPlane;
% % % % % InitialInd=find(isnan(Suite2pTable.markCycle)==1);
% % % % % 
% % % % % InitialDataI=[];
% % % % % frameGap=100;             %%%Frame before/after 1st/last whiskeCorrResults.rStim were consider to exclude to calculate neural-speed correlation.
% % % % % i=1
% % % % % tempF=double(fStim{InitialInd(i)}(:,1));
% % % % % t1=min(find(tempF>1))-frameGap;
% % % % % t2=max(find(tempF>1))+frameGap;
% % % % % ExcludeStimIndraw=[1:t1 t2:Suite2pTable.Suite2pTiffNum(InitialInd(i))/nPlane];  
% % % % % 
% % % % % ExcludeStimInd=ExcludeStimIndraw;
% % % % % for i=2:length(InitialInd)
% % % % %     ExcludeStimInd=[ExcludeStimInd,ExcludeStimInd+Suite2pTable.Suite2pTiffNum(InitialInd(i))/nPlane];
% % % % % end
% % % % % InitialInd=1:sum(Suite2pTable.Suite2pTiffNum(InitialInd))/nPlane;
% % % % % 
% % % % % 
% % % % % XTimesStdTh = 3;
% % % % % MinInterVal = 20;
% % % % % maxLag = 10;
% % % % % clear CorrResults.rSpeed pSpeed CorrResults.rStim;
% % % % % 
% % % % % for iData =1:length(NData)
% % % % % 
% % % % % for iPlane=1:nPlane
% % % % %     I1=find(CaData.CellPlaneID==iPlane);
% % % % %     [CorrResults.rSpeed(I1,1,iData),pSpeed(I1,1,iData)]=corr(NData{iData}(I1,ExcludeStimInd)',SpeedAll(ExcludeStimInd,iPlane),'type','Spearman','rows','pairwise');
% % % % % end
% % % % % 
% % % % % clear c
% % % % % StimAllT=double(StimAll>1);
% % % % % 
% % % % % 
% % % % % for iCell = 1:length(iscell)
% % % % %     iCell;
% % % % %     iPlane=CaData.CellPlaneID(iCell);
% % % % %     [c(:,iCell), lags] = xcorr(NData{iData}(iCell,InitialInd)', StimAllT(InitialInd,iPlane), maxLag, 'coeff');
% % % % %     PostI = find(lags >= 0);
% % % % %     [~, i1] = max(abs(c(PostI,iCell)));
% % % % %     CorrResults.rStim(iCell, 1,iData) = c(PostI(i1),iCell);
% % % % % end
% % % % % 
% % % % % %%Stim TTL for during anesthesia
% % % % % StimAllTAne=double(StimAllT);
% % % % % 
% % % % % 
% % % % % [CorrResults.rStimSorted,rankStimTemp]=sort(CorrResults.rStim(:, 1,1),'descend');
% % % % % figure;
% % % % % subplot('position',[0.1 0.1 0.7 0.6])
% % % % % imagesc(AmpNormalizeRow(NData{iData}(rankStimTemp,InitialInd),[1 99]))
% % % % % set(gca,'xlim',[0 length(InitialInd)],'ylim',[0 size(NData{iData},1)],'clim',[0.1 0.9])
% % % % % 
% % % % % xlabel('Time (frames)')
% % % % % ylabel('Cells')
% % % % % colormap(gray)
% % % % % b=colorbar;
% % % % % bar1=colorbar;
% % % % % bar1.Location='northoutside';
% % % % % bar1.Position=[0.4 0.92 0.2 0.01];
% % % % % bar1.Label.String=Nlabel{iData};
% % % % % % bar1.Ticks=union(ResponseLim,0);
% % % % % 
% % % % % 
% % % % % subplot('position',[0.1 0.75 0.7 0.1])
% % % % % plot(mean(StimAllT(InitialInd,:),2)>0.1)
% % % % % % bar1.Ticks=union(ResponseLim,0);
% % % % % set(gca,'xlim',[0 length(InitialInd)],'xticklabel',[],'Box','off');
% % % % % ylabel('Whisker-Stim')
% % % % % 
% % % % % 
% % % % % 
% % % % % subplot('position',[0.85 0.1 0.1 0.6])
% % % % % plot(CorrResults.rStimSorted,'color',[0.01 0.01 0.1]);
% % % % % view(90, 90);
% % % % % set(gca,'xticklabel',[],'Box','off','xlim',[0 size(NData{iData},1)]);
% % % % % ylabel('CorrResults.rStim Corr.')
% % % % % papersizePX=[0 0 16 16];
% % % % % set(gcf, 'PaperUnits', 'centimeters');
% % % % % set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
% % % % % saveas(gcf,[ResultFolder 'SponStim' Nlabel{iData}  'Corr'],'png');
% % % % % close all
% % % % % 
% % % % % 
% % % % % 
% % % % % 
% % % % % 
% % % % % [CorrResults.rSpeedSorted,rankSpeedTemp]=sort(CorrResults.rSpeed(:, 1,1),'descend');
% % % % % 
% % % % % figure;
% % % % % subplot('position',[0.1 0.1 0.7 0.6])
% % % % % imagesc(AmpNormalizeRow(NData{iData}(rankSpeedTemp,InitialInd),[1 99]))
% % % % % set(gca,'xlim',[0 length(InitialInd)],'ylim',[0 size(NData{iData},1)],'clim',[0.1 0.9])
% % % % % 
% % % % % xlabel('Time (frames)')
% % % % % ylabel('Cells')
% % % % % colormap(gray)
% % % % % b=colorbar;
% % % % % bar1=colorbar;
% % % % % bar1.Location='northoutside';
% % % % % bar1.Position=[0.4 0.92 0.2 0.01];
% % % % % bar1.Label.String=Nlabel{iData};
% % % % % 
% % % % % subplot('position',[0.1 0.75 0.7 0.1])
% % % % % plot(mean(SpeedAll(InitialInd,:),2))
% % % % % % bar1.Ticks=union(ResponseLim,0);
% % % % % set(gca,'xlim',[0 length(InitialInd)],'xticklabel',[],'Box','off');
% % % % % ylabel('Speed')
% % % % % 
% % % % % subplot('position',[0.85 0.1 0.1 0.6])
% % % % % plot(CorrResults.rSpeedSorted,'color',[0.01 0.01 0.1]);
% % % % % view(90, 90);
% % % % % set(gca,'xticklabel',[],'Box','off','xlim',[0 size(NData{iData},1)]);
% % % % % ylabel('CorrResults.rSpeed Corr.')
% % % % % papersizePX=[0 0 16 16];
% % % % % set(gcf, 'PaperUnits', 'centimeters');
% % % % % set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
% % % % % saveas(gcf,[ResultFolder 'SponSpeed' Nlabel{iData}  'Corr'],'png');
% % % % % close all
% % % % % 
% % % % % 
% % % % % % AneI=Suite2pTable.suite2pInd(Suite2pTable.AwakeState==2&Suite2pTable.PowerZero==0)
% % % % % % StimAllTAne(AneI,:)=3;
% % % % % % StimAllTAne=double(StimAllTAne==3);
% % % % % StimSLM0TAne=double(StimAllT);
% % % % % t1=Suite2pTable.suite2pInd(Suite2pTable.AwakeState==2&Suite2pTable.PowerZero==1);
% % % % % StimSLM0TAne(t1,:)=3;
% % % % % StimSLM0TAne=double(StimSLM0TAne==3);
% % % % % 
% % % % % %%Stim TTL for during awake
% % % % % StimSLM0TWake=double(StimAllT);
% % % % % t1=Suite2pTable.suite2pInd(Suite2pTable.AwakeState==1&Suite2pTable.PowerZero==1);
% % % % % StimSLM0TWake(t1,:)=3;
% % % % % StimSLM0TWake=double(StimSLM0TWake==3);
% % % % % 
% % % % % 
% % % % % 
% % % % % 
% % % % % FunExpStartI=min(find(SessFileTable.Group>0))-1;
% % % % % FunExpEndI=max(find(SessFileTable.AwakeState==1));
% % % % % FunExpInd=sum(SessFileTable.Suite2pTiffNum(1:FunExpStartI))/nPlane+1:sum(SessFileTable.Suite2pTiffNum(1:FunExpEndI))/nPlane;
% % % % % for iPlane=1:nPlane
% % % % %     I1=find(CaData.CellPlaneID==iPlane);
% % % % %     [CorrResults.rSpeed(I1,2,iData),pSpeed(I1,2,iData)]=corr(NData{iData}(I1,FunExpInd)',SpeedAll(FunExpInd,iPlane),'type','Spearman','rows','pairwise');
% % % % % end
% % % % % 
% % % % % StimAllT=StimAll>1;
% % % % % clear c;
% % % % % for iCell = 1:length(iscell)
% % % % %     iCell;
% % % % %     iPlane=CaData.CellPlaneID(iCell);
% % % % %     [c(:,iCell), lags] = xcorr(NData{iData}(iCell,FunExpInd)', StimSLM0TWake(FunExpInd,iPlane), maxLag, 'coeff');
% % % % %     PostI = find(lags >= 0);
% % % % %     [~, i1] = max(abs(c(PostI,iCell)));
% % % % %     CorrResults.rStim(iCell, 2,iData) = c(PostI(i1),iCell);
% % % % % end
% % % % % 
% % % % % 
% % % % % 
% % % % % % AneExpStartI=min(find(SessFileTable.AwakeState==2))-1;
% % % % % % AneExpEndI=max(find(SessFileTable.AwakeState==2));
% % % % % % AneExpInd=sum(SessFileTable.Suite2pTiffNum(1:AneExpStartI))/nPlane+1:sum(SessFileTable.Suite2pTiffNum(1:AneExpEndI))/nPlane;
% % % % % % 
% % % % % % clear c;
% % % % % % for iCell = 1:length(iscell)
% % % % % %     iCell;
% % % % % %     iPlane=CaData.CellPlaneID(iCell);
% % % % % %     [c(:,iCell), lags] = xcorr(NData{iData}(iCell,AneExpInd)', StimSLM0TAne(AneExpInd,iPlane), maxLag, 'coeff');
% % % % % %     PostI = find(lags >= 0);
% % % % % %     [~, i1] = max(abs(c(PostI,iCell)));
% % % % % %     CorrResults.rStim(iCell,3,iData) = c(PostI(i1),iCell);
% % % % % % end
% % % % % % 
% % % % % 
% % % % % 
% % % % % 
% % % % % 
% % % % %    P.xLeft=0.08;        %%%%%%Left Margin
% % % % %    P.xRight=0.02;       %%%%%%Right Margin
% % % % %    P.yTop=0.04;         %%%%%%Top Margin
% % % % %    P.yBottom=0.12;      %%%%%%Bottom Margin
% % % % %    P.xInt=0.08;         %%%%%%Width-interval between subplots
% % % % %    P.yInt=0.02;         %%%%%%Height-interval between subplots
% % % % % 
% % % % % 
% % % % % figure;
% % % % % for iSession = 1:size(CorrResults.rSpeed,2)
% % % % % regParam.Color=[0.1 0.1 0.1];
% % % % % regParam.Marker='o';
% % % % % regParam.MarkerSize=3;
% % % % % regParam.Rtype='pearson';
% % % % % regParam.xLim=[-0.3 0.6];
% % % % % regParam.yLim=[-0.3 0.6];
% % % % % regParam.xLabel='SpeedCorr Pre';
% % % % % regParam.yLabel='SpeedCorr Post';
% % % % % subplotLU(1,size(CorrResults.rSpeed,2),1,1,P)
% % % % % [OutPut,r,p]=LuPairRegressPlot(CorrResults.rSpeed(:,1,iData),CorrResults.rSpeed(:,2,iData),regParam);
% % % % % plot([-0.3 0.6],[-0.3 0.6],'k:');hold on;
% % % % % % set(gca,'xlim',[-0.3 0.6],'ylim',[-0.3 0.6])
% % % % % 
% % % % % regParam.xLabel='StimCorr Pre';
% % % % % regParam.yLabel='StimCorr Awake';
% % % % % regParam.xLim=[-0.05 0.2];
% % % % % regParam.yLim=[-0.05 0.2];
% % % % % subplotLU(1,size(CorrResults.rSpeed,2),1,2,P)
% % % % % [OutPut,r,p]=LuPairRegressPlot(CorrResults.rStim(:,1,iData),CorrResults.rStim(:,2,iData),regParam)    
% % % % % plot([-0.3 0.6],[-0.3 0.6],'k:');hold on;
% % % % % 
% % % % % % subplotLU(1,3,1,3,P)
% % % % % % regParam.xLabel='StimCorr Pre';
% % % % % % regParam.yLabel='StimCorr Ane';
% % % % % 
% % % % % % [OutPut,r,p]=LuPairRegressPlot(CorrResults.rStim(:,1,iData),CorrResults.rStim(:,3,iData),regParam)    
% % % % % % plot([-0.3 0.6],[-0.3 0.6],'k:');hold on;
% % % % % papersizePX=[0 0 20 10];
% % % % % set(gcf, 'PaperUnits', 'centimeters');
% % % % % set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
% % % % % saveas(gcf,[ResultFolder 'Beh' Nlabel{iData} 'CorrNoSLM'],'png');
% % % % % close all
% % % % % end
% % % % % 
% % % % % 
% % % % % 
% % % % % end
% % % % % 
% % % % % 
% % % % % figure;
% % % % % regParam.Color=[0.1 0.1 0.1];
% % % % % regParam.Marker='o';
% % % % % regParam.MarkerSize=3;
% % % % % regParam.Rtype='pearson';
% % % % % regParam.xLim=[-0.3 0.6];
% % % % % regParam.yLim=[-0.3 0.6];
% % % % % regParam.xLabel='SpeedCorr Pre DeltaF';
% % % % % regParam.yLabel='SpeedCorr Pre spks';
% % % % % subplotLU(1,2,1,1,P)
% % % % % [OutPut,r,p]=LuPairRegressPlot(CorrResults.rSpeed(:,1,1),CorrResults.rSpeed(:,1,2),regParam);
% % % % % plot([-0.3 0.6],[-0.3 0.6],'k:');hold on;
% % % % % % set(gca,'xlim',[-0.3 0.6],'ylim',[-0.3 0.6])
% % % % % 
% % % % % regParam.xLabel='StimCorr Pre DeltaF';
% % % % % regParam.yLabel='StimCorr Pre spks';
% % % % % regParam.xLim=[-0.05 0.2];
% % % % % regParam.yLim=[-0.05 0.2];
% % % % % subplotLU(1,2,1,2,P)
% % % % % [OutPut,r,p]=LuPairRegressPlot(CorrResults.rStim(:,1,1),CorrResults.rStim(:,1,2),regParam)    
% % % % % plot([-0.3 0.6],[-0.3 0.6],'k:');hold on;
% % % % % 
% % % % % 
% % % % % papersizePX=[0 0 20 10];
% % % % % set(gcf, 'PaperUnits', 'centimeters');
% % % % % set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
% % % % % saveas(gcf,[ResultFolder 'Beh' Nlabel{1} Nlabel{2} 'CorrNoSLM'],'png');
% % % % % close all
% % % % % 
% % % % % 
% % % % % responseColorMap=slanCM('seismic',64);
% % % % % SLMInfoTable=Suite2pTable(Suite2pTable.suite2pInd>0,:);
% % % % % 
% % % % % PVpower=xmlPower2PVpower(SLMTestInfo.confSet.UncagingLaserPower);
% % % % % PVpower=intersect(round(PVpower),SLMInfoTable.UncagingLaserPower);
% % % % % 
% % % % % 
% % % % % AlignedSpeed=[];
% % % % % AlignedStim=[];
% % % % % 
% % % % % for iData =1:length(NData)
% % % % % 
% % % % % tempNData=NData{iData};
% % % % % tempNData=AmpNormalizeRow(tempNData,[0 100]);
% % % % % 
% % % % % 
% % % % % 
% % % % % SLMInfoTable=Suite2pTable(Suite2pTable.suite2pInd>0,:);
% % % % % AlignedtempNData=[];
% % % % % 
% % % % % for i=1:size(SLMInfoTable,1)
% % % % %     I0=SLMInfoTable.suite2pInd(i);
% % % % %     s1=I0-PSTHparam.PreSLMCal:I0-1;
% % % % %     s2=I0:I0+PSTHparam.PostSLMCal-1;
% % % % %     AlignedtempNData(:,:,i)=tempNData(:,[s1 s2])-repmat(nanmean(tempNData(:,s1),2),1,length(s1)+length(s2));
% % % % % 
% % % % %     if iData==1
% % % % %     AlignedSpeed(:,i)=mean(SpeedAll([s1 s2],:),2);
% % % % %     AlignedStim(:,i)=mean(StimAll([s1 s2],:),2);
% % % % %     end
% % % % % end
% % % % % 
% % % % % TrialThNum=3;
% % % % % 
% % % % % PSTHparam.TestStepFrame=3;
% % % % % [AlignedtempNData,TargetInfoTable,CellResponse,statCellRes,TargetResponse,TargetCellResP,TargetCellResR,CellSampleN]=Aligned_FromSuite2p(NData{iData},TargetCellList,iscell,Suite2pTable,PVpower,PSTHparam);
% % % % % if iData==1
% % % % % 
% % % % % FDR=0.1;
% % % % % [h, crit_p, adj_p]=fdr_bh(TargetCellResP((~isnan(TargetCellResP))&CellSampleN>=TrialThNum),FDR,'pdep');
% % % % % crit_p=min([crit_p 0.05]);
% % % % % 
% % % % % 
% % % % % PowerTargetI=zeros(length(TargetCellList),1);
% % % % % SuccTarget=PowerTargetI;
% % % % % for iCell=1:length(TargetCellList)
% % % % %     for iPower=1:length(PVpower)
% % % % %         if TargetCellResP(iCell,iPower)<=crit_p&&TargetCellResR(iCell,iPower)>0&&CellSampleN(iCell,iPower)>=TrialThNum
% % % % %            PowerTargetI(iCell)=iPower;
% % % % %            SuccTarget(iCell)=1;
% % % % %            SuccAmp(iCell)=TargetCellResR(iCell,iPower);
% % % % %         end
% % % % %     end
% % % % % end
% % % % % sum(SuccTarget)
% % % % % 
% % % % % pAll=[];
% % % % % for iCell=1:length(TargetCellList)
% % % % %     if SuccTarget(iCell)
% % % % %        temp=statCellRes(iCell,PowerTargetI(iCell)).p(:);
% % % % %        % temp(TargetCellList(iCell))=[];
% % % % %        pAll=[pAll;temp];
% % % % %     end
% % % % % end
% % % % % FDR=0.1;
% % % % % [h, crit_pAll, adj_p]=fdr_bh(pAll,FDR,'pdep');
% % % % % crit_pAll=min([crit_pAll 0.05]);
% % % % % 
% % % % % end
% % % % % 
% % % % % 
% % % % % 
% % % % % close all
% % % % % 
% % % % % Param.PlotType=3
% % % % % Param.statisP=0;
% % % % % Param.LegendShow=0;
% % % % % Param.Legend=[]
% % % % % 
% % % % % 
% % % % %    P.xLeft=0.06;        %%%%%%Left Margin
% % % % %    P.xRight=0.02;       %%%%%%Right Margin
% % % % %    P.yTop=0.02;         %%%%%%Top Margin
% % % % %    P.yBottom=0.1;      %%%%%%Bottom Margin
% % % % %    P.xInt=0.02;         %%%%%%Width-interval between subplots
% % % % %    P.yInt=0.01;         %%%%%%Height-interval between subplots
% % % % % 
% % % % % figure;
% % % % % for iCell = 1:length(TargetCellList)
% % % % %     for iPower=1:length(PVpower)
% % % % %         if ~isempty(TargetResponse(iCell,iPower))&CellSampleN(iCell,iPower)>=TrialThNum
% % % % % 
% % % % %             % TargetCellList(iCell)
% % % % %             subplotLU(length(TargetCellList),length(PVpower),iCell,iPower,P)
% % % % %             RateHist_GroupPlot(TimBinFrame+0.5,TargetResponse(iCell,iPower),[0.1 0.1 0.1],Param)
% % % % %             hold on;
% % % % %             text(-10,0.1,['C' num2str(iCell) 'n = ' num2str(CellSampleN(iCell,iPower))])
% % % % %             set(gca,'xlim',[-PSTHparam.PreSLMCal PSTHparam.PostSLMCal],'xtick',[-PSTHparam.PreSLMCal:5:PSTHparam.PostSLMCal])
% % % % %             ylabel(['Target' num2str(iCell)]);
% % % % %         end
% % % % %     end
% % % % % 
% % % % % end
% % % % % papersizePX=[0 0 length(PVpower)*4 5*length(TargetCellList) ];
% % % % % set(gcf, 'PaperUnits', 'centimeters');
% % % % % set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
% % % % % saveas(gcf,[ResultFolder 'AllSLMTargetResponse' Nlabel{iData}],'png');
% % % % % close all
% % % % % 
% % % % % 
% % % % % 
% % % % % figure;
% % % % % for iCell = 1:length(TargetCellList)
% % % % %     for iPower=1:length(PVpower)
% % % % %         if ~isempty(CellResponse{iCell,iPower})&CellSampleN(iCell,iPower)>=TrialThNum
% % % % %             % TargetCellList(iCell)
% % % % %             subplotLU(length(TargetCellList),length(PVpower),iCell,iPower,P)
% % % % %             imagesc(TimBinFrame+0.5,1:length(iscell),CellResponse{iCell,iPower});
% % % % %             colormap(responseColorMap);
% % % % %             set(gca,'clim',[-0.1 0.1],'ylim',[0 length(iscell)+1]);
% % % % %             hold on;
% % % % %             plot(TimBinFrame(1),TargetCellList(iCell)+0.5,'g>');
% % % % %             set(gca,'xlim',[-PSTHparam.PreSLMCal PSTHparam.PostSLMCal],'xtick',[-PSTHparam.PreSLMCal:5:PSTHparam.PostSLMCal])
% % % % %         end
% % % % %     end
% % % % % 
% % % % % end
% % % % % papersizePX=[0 0 length(PVpower)*3 5*length(TargetCellList) ];
% % % % % set(gcf, 'PaperUnits', 'centimeters');
% % % % % set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
% % % % % saveas(gcf,[ResultFolder 'AllSLMTargetResponseMap' Nlabel{iData}],'png');
% % % % % close all
% % % % % 
% % % % % 
% % % % % 
% % % % % 
% % % % % AlignedNData{iData}=AlignedtempNData;
% % % % % 
% % % % % 
% % % % % 
% % % % % 
% % % % % end
% % % % % 
% % % % % 
% % % % % save([ResultFolder 'Step1Meta.mat']);
% % % % % 
% % % % % for iData=1:length(NData)
% % % % % clear Data1 Data2
% % % % % for iFun=1:length(SLMPosInfo.Group)
% % % % %     Data1{iFun}=CorrResults.rSpeed(TargetCellList(TargetCellListFunGroup==iFun),1,iData);
% % % % %     Data1{iFun+length(SLMPosInfo.Group)}=CorrResults.rSpeed(TargetCellList(TargetCellListFunGroup==iFun),2,iData);
% % % % % 
% % % % %     Data2{iFun}=CorrResults.rStim(TargetCellList(TargetCellListFunGroup==iFun),1,iData);
% % % % %     Data2{iFun+length(SLMPosInfo.Group)}=CorrResults.rStim(TargetCellList(TargetCellListFunGroup==iFun),2,iData);
% % % % %     % Data2{iFun+length(SLMPosInfo.Group)*2}=CorrResults.rStim(TargetCellList(TargetCellListFunGroup==iFun),3,iData);
% % % % % 
% % % % % end
% % % % % 
% % % % % x1=1:length(SLMPosInfo.Group);x1=[x1 x1+length(SLMPosInfo.Group)+1];
% % % % % % x2=1:length(SLMPosInfo.Group);x2=[x2 x2+length(SLMPosInfo.Group)+1 x2+length(SLMPosInfo.Group)*2+2];
% % % % % x2=1:length(SLMPosInfo.Group);x2=[x2 x2+length(SLMPosInfo.Group)+1];
% % % % % 
% % % % % GroupLabel={'L.','S.','N.'};
% % % % % 
% % % % % nGroup=length(SLMPosInfo.Group);
% % % % % 
% % % % % 
% % % % % 
% % % % % % % subplotLU(2,1,1,1)
% % % % % % % stats=ErrorBoxPlotLU(x1,Data1,repmat(GroupColor,2,1),[ResultFolder 'SpeedGroup']);
% % % % % % % LuLegend([9 9 9;0.5 0.4 0.3;10 10 10;0.5 0.4 0.3],0,GroupLabel,GroupColor,8);
% % % % % % % set(gca,'xlim',[0 12],'ylim',[-0.2 1]);
% % % % % % % subplotLU(2,1,2,1)
% % % % % % % stats=ErrorBoxPlotLU(x2,Data2,repmat(GroupColor,3,1),[ResultFolder 'SpeedGroup']);
% % % % % % % set(gca,'xlim',[0 12]);
% % % % % 
% % % % % GroupColor=[247 150 111;239 109 249;121 247 111]/255;
% % % % % 
% % % % %    GroupPair.CorrName='fdr';
% % % % %    GroupPair.Q=0.1;
% % % % % 
% % % % % 
% % % % %    GroupPair.SignY=1;
% % % % %    GroupPair.Plot=1;
% % % % %    GroupPair.Std=1;      %%%%%%%%%using standard deviation as errorbar
% % % % %    GroupPair.SamplePlot=1; %%%%%%%%%Plot Individual Sample Point
% % % % %    GroupPair.SamplePairedPlot=1; %%%%%%%%%Dash line for paired comparison sample
% % % % %    GroupPair.LimY=[0 GroupPair.SignY*1.2];
% % % % %    GroupPair.Marker={'o'};
% % % % %    GroupPair.GroupName=repmat(GroupLabel,1,2);
% % % % %    % GroupID=repmat(1:nGroup,1,2);
% % % % % 
% % % % % P.xInt=0.1;
% % % % % P.xLeft=0.1;
% % % % % 
% % % % % figure;
% % % % % subplotLU(1,2,1,1,P)
% % % % %    p1=[1 1 2;2 3 3];p1=[p1 p1+nGroup];
% % % % %    p2=[1:nGroup];p2=[p2;p2+nGroup];
% % % % %    GroupPair.Pair=[p1 p2];
% % % % % GroupPair.GroupName=repmat(GroupLabel,1,2);
% % % % %    GroupPair.SignY=0.5;
% % % % % 
% % % % % ErrorBarPlotLU(x1,Data1,[],repmat(GroupColor,2,1),2,1,[ResultFolder Nlabel{iData} 'SpeedGroup.txt'],GroupPair,repmat(1:nGroup,1,2));
% % % % % LuLegend([9 9 9;0.5 0.4 0.3;10 10 10;0.5 0.4 0.3],0,GroupLabel,GroupColor,8);
% % % % % set(gca,'xlim',[0 12],'ylim',[-0.2 0.6],'ytick',[-0.2:0.2:0.6],'xtick',[2 6],'xticklabel',{'Spon 1','Spon 2'});
% % % % % ylabel('Speed Corr')
% % % % % 
% % % % % subplotLU(1,2,1,2,P)
% % % % %    % p1=[1 1 2;2 3 3];p1=[p1 p1+nGroup p1+nGroup*2];
% % % % %    % p2=[1:nGroup];p2=[p2 p2 p2+nGroup;p2+nGroup p2+nGroup*2 p2+nGroup*2];
% % % % %   p1=[1 1 2;2 3 3];p1=[p1 p1+nGroup];
% % % % %    p2=[1:nGroup];p2=[p2;p2+nGroup];
% % % % % 
% % % % %    GroupPair.Pair=[p1 p2];
% % % % %    GroupPair.SignY=0.25;
% % % % %    GroupID=repmat(1:nGroup,1,2);
% % % % %    GroupPair.GroupName=repmat(GroupLabel,1,2);
% % % % % 
% % % % % ErrorBarPlotLU(x2,Data2,[],repmat(GroupColor,2,1),2,1,[ResultFolder Nlabel{iData} 'StimGroup.txt'],GroupPair,repmat(1:nGroup,1,2));
% % % % % set(gca,'xlim',[0 12],'ylim',[-0.05 0.3],'ytick',[-0.05:0.05:0.3],'xtick',[2 6],'xticklabel',{'Spon 1','Spon 2'});
% % % % % ylabel('Stim Corr')
% % % % % 
% % % % % papersizePX=[0 0 22 10];
% % % % % set(gcf, 'PaperUnits', 'centimeters');
% % % % % set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
% % % % % saveas(gcf,[ResultFolder 'BehCorrCompare' Nlabel{iData}],'png');
% % % % % close all
% % % % % 
% % % % % end
% % % % % 
% % % % % 
% % % % % end
