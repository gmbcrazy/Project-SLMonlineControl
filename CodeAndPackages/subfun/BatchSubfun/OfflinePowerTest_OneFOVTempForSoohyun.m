function OfflinePowerTest_OneFOVTempForSoohyun(FOVtemp, Suite2pDataKeywords, suite2pFOVPathLocalTemp,PSTHparam,SaveFunCon)


% GroupLabel={'L','S','N'};
% nGroup=length(GroupLabel);
% GroupColor=[255 51 153;91 20 212;121 247 111]/255;
% NodeColor=repmat([0.9 0.9 0.9],length(iscell),1);
% 
papersizePX=[0 0 33 9];
papersizePXX=[0 0 22 16];

% suite2pFOVPathLocalTemp=suite2pFOVPathLocal{iFOV};
resultPaths = findAllFoldersKeyWords(suite2pFOVPathLocalTemp, Suite2pDataKeywords);

metaDataPath=Get_ExpDataFolder(resultPaths{1},'Step1Basic',{'Step1Meta.mat','BehAll.mat','.png'});
load([metaDataPath 'Step1Meta.mat'])
confSet=SLMPosInfo.confSetFinal
SLMPos3D=SLMTestInfo.Pos3Dneed;
TimBinFrame=-PSTHparam.PreSLMCal:PSTHparam.PostSLMCal-1;        %<----------------------------------------------------------------------------------Edit,Frame # before SLM to calculate baseline map

Zdepth = confSet.scan_Z + confSet.ETL;
suite2pPath = findAllFoldersKeyWords(suite2pFOVPathLocalTemp, Suite2pDataKeywords);

confSet.save_path0=suite2pPath{1};
[NeuronPos3D,NeuronPos3DRaw,CaData,CaDataPlane,Neuronstat]=Extract_Suite2p(confSet);
% load('C:\Users\zhangl33\Projects\Project-SLMonlineControl\CodeAndPackages\subfun\Color\colorMapPN3.mat','colorMapPN1')
[~, ~, CellMedCenter, cellBoundary, ~] = Suite2pCellIDMapFromStat(CaData.statCell, [confSet.SLM_Pixels_Y confSet.SLM_Pixels_X]);
Cell3DPos=[CellMedCenter Zdepth(CaData.CellPlaneID)'];
CellDist=squareform(pdist(Cell3DPos));

iData=PSTHparam.iData;
TestStepFrame=PSTHparam.TestStepFrame;
iscell=find(CaData.iscell(CaData.iscell(:,1)>0));
PVpower=xmlPower2PVpower(SLMTestInfo.confSet.UncagingLaserPower);
PVpower=intersect(round(PVpower),SLMInfoTable.UncagingLaserPower);
GroupLabel={'L','S','N'};
nGroup=length(GroupLabel);
GroupColor=[255 51 153;91 20 212;121 247 111]/255;
NodeColor=repmat([0.9 0.9 0.9],length(iscell),1);



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
ParamNet.xMat=[0.01 0.15 0.15 0.15 0.15 0.15];
ParamNet.yMat=[0.7 0.7 0.7 0.7 0.7 0.7];


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
       ResonseSpeedDep=[ResonseSpeedDep;squeeze(CorrResults.rSpeed(:,1,iData))];
       ResonseStimDep=[ResonseStimDep;squeeze(CorrResults.rStim(:,1,iData))];


       DistTargetDep = [DistTargetDep; CellDist(TargetCellList(iCell),:)'];
       
       tempDist=CellDist(TargetCellList(iCell),:);
       NonSamePlaneI=find(abs(NeuronPos3D(:,3)-NeuronPos3D(TargetCellList(iCell),3)))>0.1;
       tempDist(NonSamePlaneI)=NaN;
       DistPlaneTargetDep=[DistPlaneTargetDep;tempDist(:)];
    end
end


MixedData =[DataResponse DistTargetDep DistPlaneTargetDep StimFunGDep ResonseSpeedDep ResonseStimDep TargetNameDep];
% MixedData(:,end+1) = iFOV;

tbl = table(DataResponse, DistTargetDep, DistPlaneTargetDep, StimFunGDep, ...
    ResonseSpeedDep, ResonseStimDep, TargetNameDep,  ...
    'VariableNames', {'Response', 'Dist','DistPlane','StimGroup','SpeedR','StimR','Cell'});


% temptbl=tbl;
% temptbl(Invalid,:)=[];


[~,Session,~]=fileparts(suite2pFOVPathLocalTemp(1:end-1));

SaveP1=[SaveFunCon Session '\'];
mkdir(SaveP1);
EGTestFolder=[SaveP1 'ExampleSLMTest\'];
mkdir(EGTestFolder)


% csvwrite([SaveP1 'tempMixedLinearM.csv'], 'tbl');
writetable(tbl,[SaveFunCon 'PowerTestResponse.csv']);


DistScoreTemp=sum(DistTarget)/sum(SuccTarget)*SLMPosInfo.yaml.umPerlPixelX;

ParamNetDist=ParamNet;
ParamNetDist.AllCellRes=1; 
ParamNetDist.ScoreLim=[0;450];
ParamNetDist.ScoreMap=slanCM('amethyst',64);

PowerTestAdj=PowerTestAdj-diag(diag(PowerTestAdj));


close all

% figure;
ParamNetDist.AllCellRes=1; 
ParamNetDist.ScoreLabel='Dist. from Target (µm)';
[GgraphOut,Res,r,p]=ResPosNetwork_ScoreNodeOrder(SLMPosInfo,PowerTestAdj,PowerTestNode,ParamNetDist,DistScoreTemp(:))
% GgraphOut.G.Nodes.Weight(GgraphOut.G.Nodes.Weight>0)=10;
GgraphOut.p.ArrowSize=4;
GgraphOut.p.MarkerSize(GgraphOut.G.Nodes.Weight>0)=2;

% GgraphOut.p.MarkerSize=2;
GgraphOut.p.LineWidth=0.5;
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
LuFontStandard
saveas(gcf,[SaveP1 'SingleSLMVsDist_SoohyunPresnetation'],'tif');
% print(gcf, [SaveP1 'SingleSLMVsDist.svg'], '-dsvg', '-painters');


ParamNetDist.AllCellRes=0; 
for iFun=1:length(GroupLabel)

    tempAdj=PowerTestAdj;
    tempAdj(TargetCellList(TargetCellListFunGroup~=iFun),:)=NaN;

    tempNode=PowerTestNode;
    tempNode(TargetCellList(TargetCellListFunGroup~=iFun))=0;
    [GgraphOut,Res,r,p]=ResPosNetwork_ScoreNodeOrder(SLMPosInfo,tempAdj,tempNode,ParamNetDist,DistScoreTemp(:))
GgraphOut.p.ArrowSize=4;
GgraphOut.p.MarkerSize(GgraphOut.G.Nodes.Weight>0)=2;
GgraphOut.p.LineWidth=0.5;

set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
LuFontStandard
% print(gcf, [SaveP1 GroupLabel{iFun} 'SingleSLMVsDist.svg'], '-dsvg', '-painters');
% print(gcf, [SaveP1 GroupLabel{iFun} 'SingleSLMVsSpeed.eps'], '-depsc', '-painters');
print(gcf, [SaveP1 GroupLabel{iFun} 'SingleSLMVsDist_SoohyunPresnetation.tif'], '-dtiffn', '-painters');
% print(gcf, [SaveP1 'SingleSLMVsSpeed.png'], '-dpng', '-painters');
end


close all







