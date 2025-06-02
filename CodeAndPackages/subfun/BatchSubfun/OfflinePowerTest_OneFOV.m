function OfflinePowerTest_OneFOV(FOVtemp, Suite2pDataKeywords, suite2pFOVPathLocalTemp,PSTHparam,SaveFunCon)


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

Zdepth = confSet.scan_Z + confSet.ETL;
suite2pPath = findAllFoldersKeyWords(suite2pFOVPathLocalTemp, Suite2pDataKeywords);

confSet.save_path0=suite2pPath{1};
[NeuronPos3D,NeuronPos3DRaw,CaData,CaDataPlane,Neuronstat]=Extract_Suite2p(confSet);
% load('C:\Users\zhangl33\Projects\Project-SLMonlineControl\CodeAndPackages\subfun\Color\colorMapPN3.mat','colorMapPN1')
[~, ~, CellMedCenter, cellBoundary, ~] = Suite2pCellIDMapFromStat(CaData.statCell, [confSet.SLM_Pixels_Y confSet.SLM_Pixels_X]);
Cell3DPos=[CellMedCenter Zdepth(CaData.CellPlaneID)'];
CellDist=squareform(pdist(Cell3DPos));

iData=PSTHparam.iData;
PSTHparam.TestStepFrame=3;
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

% csvwrite([SaveP1 'tempMixedLinearM.csv'], 'tbl');
writetable(tbl,[SaveFunCon 'PowerTestResponse.csv']);


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
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
LuFontStandard
saveas(gcf,[SaveP1 'SingleSLMVsDist'],'tif');


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

set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
LuFontStandard
% print(gcf, [SaveP1 GroupLabel{iFun} 'SingleSLMVsDist.svg'], '-dsvg', '-painters');
% print(gcf, [SaveP1 GroupLabel{iFun} 'SingleSLMVsSpeed.eps'], '-depsc', '-painters');
print(gcf, [SaveP1 GroupLabel{iFun} 'SingleSLMVsDist.tif'], '-dtiffn', '-painters');
% print(gcf, [SaveP1 'SingleSLMVsSpeed.png'], '-dpng', '-painters');
end







rSpeed=CorrResults.rSpeed;
rStim=CorrResults.rStim;
ParamNet.AllCellRes=1; 
close all
[GgraphOut,Res,r,p]=ResPosNetwork_ScoreNodeOrder(SLMPosInfo,PowerTestAdj,PowerTestNode,ParamNet,rSpeed(:,1,iData))
GgraphOut.p.ArrowSize=4;
% GgraphOut.p.MarkerSize=2;
GgraphOut.p.LineWidth=0.5;
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
% saveas(gcf,[SaveP1 'SingleSLMVsSpeed'],'tif');
% saveas(gcf,[SaveP1 'SingleSLMVsSpeed.eps'],'epsc');
% saveas(gcf,[SaveP1 'SingleSLMVsSpeed'],'fig');
LuFontStandard
% print(gcf, [SaveP1 'SingleSLMVsSpeed.svg'], '-dsvg', '-painters');
% print(gcf, [SaveP1 'SingleSLMVsSpeed.eps'], '-depsc', '-painters');
print(gcf, [SaveP1 'SingleSLMVsSpeed.tif'], '-dtiffn', '-painters');
% print(gcf, [SaveP1 'SingleSLMVsSpeed.png'], '-dpng', '-painters');
close all
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

set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
LuFontStandard
% print(gcf, [SaveP1 GroupLabel{iFun} 'SingleSLMVsSpeed.svg'], '-dsvg', '-painters');
% print(gcf, [SaveP1 GroupLabel{iFun} 'SingleSLMVsSpeed.eps'], '-depsc', '-painters');
print(gcf, [SaveP1 GroupLabel{iFun} 'SingleSLMVsSpeed.tif'], '-dtiffn', '-painters');
% print(gcf, [SaveP1 'SingleSLMVsSpeed.png'], '-dpng', '-painters');
% 
% 
close all

end






ParamNet.AllCellRes=0; 
P.xLeft=0.1;
P.yBottom=0.12;
P.yTop=0.04;
P.xInt=0.04;
P.yInt=0.1;
P.xRight=0.1;

ParamNet.AllCellRes=1; 

ParamNet.ScoreLim=[-0.1 0.1];
ParamNet.ScoreLabel='Stim Corr.';
[GgraphOut,Res,r,p]=ResPosNetwork_ScoreNodeOrder(SLMPosInfo,PowerTestAdj,PowerTestNode,ParamNet,rStim(:,1,iData))
GgraphOut.p.ArrowSize=8;
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
% print(gcf, [SaveP1 'SingleSLMVsStim.svg'], '-dsvg', '-painters');
% print(gcf, [SaveP1 'SingleSLMVsStim.eps'], '-depsc', '-painters');
print(gcf, [SaveP1 'SingleSLMVsStim.tif'], '-dtiff', '-painters');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
LuFontStandard
% print(gcf, [SaveP1 'SepFunSingleSLMVsSpeedRegress.svg'], '-dsvg', '-painters');
% print(gcf, [SaveP1 'SepFunSingleSLMVsSpeedRegress.eps'], '-depsc', '-painters');
print(gcf, [SaveP1 'SepFunSingleSLMVsStimRegress.tif'], '-dtiffn', '-painters');
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

set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
LuFontStandard
% print(gcf, [SaveP1 GroupLabel{iFun} 'SingleSLMVsStim.svg'], '-dsvg', '-painters');
% print(gcf, [SaveP1 GroupLabel{iFun} 'SingleSLMVsSpeed.eps'], '-depsc', '-painters');
print(gcf, [SaveP1 GroupLabel{iFun} 'SingleSLMVStim.tif'], '-dtiffn', '-painters');
% print(gcf, [SaveP1 'SingleSLMVsSpeed.png'], '-dpng', '-painters');
% 
% 
close all

end



P.xLeft=0.1;
P.yBottom=0.12;
P.yTop=0.04;
P.xInt=0.04;
P.yInt=0.1;
ParamNetDist.AllCellRes=0;
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
   ResParam.xLim=ParamNetDist.ScoreLim;
   ResParam.yLim=[0 ParamNetDist.ResponseLim(2)];
   ResParam.xLabel=ParamNetDist.ScoreLabel;

   if iFun==1
   ResParam.yLabel='Excitation';
   else
   ResParam.yLabel='';
   end


subplotLU(2,length(GroupLabel),1,iFun,P);
if ParamNetDist.AllCellRes
[Res,r,p]=LuPairRegressPlot(squeeze(DistScoreTemp)',PosRate(:),ResParam);
else
[Res,r,p]=LuPairRegressPlot(squeeze(DistScoreTemp(PosRate>0))',PosRate(PosRate>0),ResParam);
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
ResParam.yLim=[ParamNetDist.ResponseLim(1) 0];

subplotLU(2,length(GroupLabel),2,iFun,P);
if ParamNetDist.AllCellRes
[Res,r,p]=LuPairRegressPlot(squeeze(DistScoreTemp)',NegRate(:),ResParam);
else
[Res,r,p]=LuPairRegressPlot(squeeze(DistScoreTemp(NegRate<0))',NegRate(NegRate<0),ResParam);
end
   if iFun~=1
   set(gca,'yticklabel',[])
   end
end
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePXX,'PaperSize',papersizePXX(3:4));
LuFontStandard
print(gcf, [SaveP1 'SepFunSingleSLMVsDistRegress.tif'], '-dtiffn', '-painters');
% print(gcf, [ResultFolder 'SingleSLMVsSpeed.png'], '-dpng', '-painters');



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
[Res,r,p]=LuPairRegressPlot(squeeze(rSpeed(:,1,2)),PosRate(:),ResParam);
else
[Res,r,p]=LuPairRegressPlot(squeeze(rSpeed(PosRate>0,1,2)),PosRate(PosRate>0),ResParam);
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
[Res,r,p]=LuPairRegressPlot(squeeze(rSpeed(:,1,2)),NegRate(:),ResParam);
else
[Res,r,p]=LuPairRegressPlot(squeeze(rSpeed(NegRate<0,1,2)),NegRate(NegRate<0),ResParam);
end
   if iFun~=1
   set(gca,'yticklabel',[])
   end
end
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePXX,'PaperSize',papersizePXX(3:4));
LuFontStandard
print(gcf, [SaveP1 'SepFunSingleSLMVsSpeedRegress.tif'], '-dtiffn', '-painters');
% print(gcf, [ResultFolder 'SingleSLMVsSpeed.png'], '-dpng', '-painters');



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
% papersizePX=[0 0 22 16];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePXX,'PaperSize',papersizePXX(3:4));
LuFontStandard
% print(gcf, [SaveP1 'SepFunSingleSLMVsStimRegress.svg'], '-dsvg', '-painters');
% print(gcf, [SaveP1 'SepFunSingleSLMVsStimRegress.eps'], '-depsc', '-painters');
print(gcf, [SaveP1 'SepFunSingleSLMVsStimRegress.tif'], '-dtiffn', '-painters');




