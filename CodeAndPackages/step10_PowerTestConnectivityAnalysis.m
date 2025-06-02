clear all
WorkFolder='E:\LuSLMOnlineTest\SL0855-Emx1G6CII-AAV9CAMKII\03042025\';
ResultFolder = Get_ExpDataFolder(WorkFolder, 'Step1', {'Step1Meta.mat','.png'})
load([ResultFolder 'Step1Meta.mat'])
% load('C:\Users\zhangl33\Projects\Project-SLMonlineControl\CodeAndPackages\subfun\Color\colorMapPN3.mat','colorMapPN1')
[~, ~, CellMedCenter, cellBoundary, ~] = Suite2pCellIDMapFromStat(CaData.statCell, [confSet.SLM_Pixels_Y confSet.SLM_Pixels_X]);
Cell3DPos=[CellMedCenter Zdepth(CaData.CellPlaneID)'];
CellDist=squareform(pdist(Cell3DPos));



iData=1;
PSTHparam.TestStepFrame=3;
% [AlignedtempNData,AlignedInfoTable,statCellRes,TargetResponse,TargetCellResP,TargetCellResR,CellSampleN]=Aligned_FromSuite2p(NData{iData},TargetCellList,iscell,Suite2pTable,PVpower,PSTHparam);
[AlignedtempNData,AlignedInfoTable,CellResponse,statCellRes,TargetResponse,TargetCellResP,TargetCellResR,CellSampleN]=Aligned_FromSuite2p(NData{iData},TargetCellList,iscell,Suite2pTable,PVpower,PSTHparam);


ResultFolder=[ProcessFolder 'Results\'];
mkdir(ResultFolder)
ResultFolder=[ResultFolder 'Step2\'];


EGTestFolder=[ResultFolder 'ExampleSLMTest\'];
mkdir(EGTestFolder)

mkdir(ResultFolder)
SubPlotsFolder=[ResultFolder 'SubPlots\'];
mkdir(SubPlotsFolder)


GroupLabel={'L','S','N'};
nGroup=length(SLMPosInfo.Group);
GroupColor=[255 51 153;91 20 212;121 247 111]/255;


% TargetCellBand=zeros(length(iscell),1);
% TargetCellBand(TargetCellList)=TargetCellListFunGroup;
% 
NodeColor=repmat([0.9 0.9 0.9],length(iscell),1);

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
% saveas(gcf,[ResultFolder 'SingleSLMVsSpeed'],'tif');
% saveas(gcf,[ResultFolder 'SingleSLMVsSpeed.eps'],'epsc');
% saveas(gcf,[ResultFolder 'SingleSLMVsSpeed'],'fig');
LuFontStandard
print(gcf, [ResultFolder 'SingleSLMVsDist.svg'], '-dsvg', '-painters');
% print(gcf, [ResultFolder 'SingleSLMVsDist.eps'], '-depsc', '-painters');
print(gcf, [ResultFolder 'SingleSLMVsDist.tif'], '-dtiffn', '-painters');
% print(gcf, [ResultFolder 'SingleSLMVsSpeed.png'], '-dpng', '-painters');

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
print(gcf, [ResultFolder GroupLabel{iFun} 'SingleSLMVsDist.svg'], '-dsvg', '-painters');
% print(gcf, [ResultFolder GroupLabel{iFun} 'SingleSLMVsSpeed.eps'], '-depsc', '-painters');
print(gcf, [ResultFolder GroupLabel{iFun} 'SingleSLMVsDist.tif'], '-dtiffn', '-painters');
% print(gcf, [ResultFolder 'SingleSLMVsSpeed.png'], '-dpng', '-painters');
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
papersizePX=[0 0 22 16];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
LuFontStandard
print(gcf, [ResultFolder 'SepFunSingleSLMVsSpeedRegress.svg'], '-dsvg', '-painters');
print(gcf, [ResultFolder 'SepFunSingleSLMVsSpeedRegress.eps'], '-depsc', '-painters');
print(gcf, [ResultFolder 'SepFunSingleSLMVsSpeedRegress.tif'], '-dtiffn', '-painters');
% print(gcf, [ResultFolder 'SingleSLMVsSpeed.png'], '-dpng', '-painters');











ParamNet.AllCellRes=1; 
close all
[GgraphOut,Res,r,p]=ResPosNetwork_ScoreNodeOrder(SLMPosInfo,PowerTestAdj,PowerTestNode,ParamNet,rSpeed(:,1,2))
GgraphOut.p.ArrowSize=4;
% GgraphOut.p.MarkerSize=2;
GgraphOut.p.LineWidth=0.5;
papersizePX=[0 0 33 9];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
% saveas(gcf,[ResultFolder 'SingleSLMVsSpeed'],'tif');
% saveas(gcf,[ResultFolder 'SingleSLMVsSpeed.eps'],'epsc');
% saveas(gcf,[ResultFolder 'SingleSLMVsSpeed'],'fig');
LuFontStandard
print(gcf, [ResultFolder 'SingleSLMVsSpeed.svg'], '-dsvg', '-painters');
print(gcf, [ResultFolder 'SingleSLMVsSpeed.eps'], '-depsc', '-painters');
print(gcf, [ResultFolder 'SingleSLMVsSpeed.tif'], '-dtiffn', '-painters');
% print(gcf, [ResultFolder 'SingleSLMVsSpeed.png'], '-dpng', '-painters');

ParamNet.AllCellRes=0; 
for iFun=1:length(GroupLabel)

    tempAdj=PowerTestAdj;
    tempAdj(TargetCellList(TargetCellListFunGroup~=iFun),:)=NaN;

    tempNode=PowerTestNode;
    tempNode(TargetCellList(TargetCellListFunGroup~=iFun))=0;
    [GgraphOut,Res,r,p]=ResPosNetwork_ScoreNodeOrder(SLMPosInfo,tempAdj,tempNode,ParamNet,rSpeed(:,1,2))
GgraphOut.p.ArrowSize=4;
% GgraphOut.p.MarkerSize=2;
GgraphOut.p.LineWidth=0.5;

papersizePX=[0 0 33 9];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
LuFontStandard
print(gcf, [ResultFolder GroupLabel{iFun} 'SingleSLMVsSpeed.svg'], '-dsvg', '-painters');
% print(gcf, [ResultFolder GroupLabel{iFun} 'SingleSLMVsSpeed.eps'], '-depsc', '-painters');
print(gcf, [ResultFolder GroupLabel{iFun} 'SingleSLMVsSpeed.tif'], '-dtiffn', '-painters');
% print(gcf, [ResultFolder 'SingleSLMVsSpeed.png'], '-dpng', '-painters');
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
papersizePX=[0 0 22 16];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
LuFontStandard
print(gcf, [ResultFolder 'SepFunSingleSLMVsSpeedRegress.svg'], '-dsvg', '-painters');
print(gcf, [ResultFolder 'SepFunSingleSLMVsSpeedRegress.eps'], '-depsc', '-painters');
print(gcf, [ResultFolder 'SepFunSingleSLMVsSpeedRegress.tif'], '-dtiffn', '-painters');
% print(gcf, [ResultFolder 'SingleSLMVsSpeed.png'], '-dpng', '-painters');




ParamNet.AllCellRes=1; 

ParamNet.ScoreLim=[-0.1 0.1];
ParamNet.ScoreLabel='Stim Corr.';
[GgraphOut,Res,r,p]=ResPosNetwork_ScoreNodeOrder(SLMPosInfo,PowerTestAdj,PowerTestNode,ParamNet,rStim(:,1,2))
GgraphOut.p.ArrowSize=8;
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
print(gcf, [ResultFolder 'SingleSLMVsStim.svg'], '-dsvg', '-painters');
print(gcf, [ResultFolder 'SingleSLMVsStim.eps'], '-depsc', '-painters');
print(gcf, [ResultFolder 'SingleSLMVsStim.tif'], '-dtiff', '-painters');


ParamNet.AllCellRes=0; 
for iFun=1:length(GroupLabel)

    tempAdj=PowerTestAdj;
    tempAdj(TargetCellList(TargetCellListFunGroup~=iFun),:)=NaN;

    tempNode=PowerTestNode;
    tempNode(TargetCellList(TargetCellListFunGroup~=iFun))=0;
    [GgraphOut,Res,r,p]=ResPosNetwork_ScoreNodeOrder(SLMPosInfo,tempAdj,tempNode,ParamNet,rStim(:,1,2))
GgraphOut.p.ArrowSize=4;
% GgraphOut.p.MarkerSize=2;
GgraphOut.p.LineWidth=0.5;

papersizePX=[0 0 33 9];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
LuFontStandard
print(gcf, [ResultFolder GroupLabel{iFun} 'SingleSLMVsStim.svg'], '-dsvg', '-painters');
% print(gcf, [ResultFolder GroupLabel{iFun} 'SingleSLMVsSpeed.eps'], '-depsc', '-painters');
print(gcf, [ResultFolder GroupLabel{iFun} 'SingleSLMVStim.tif'], '-dtiffn', '-painters');
% print(gcf, [ResultFolder 'SingleSLMVsSpeed.png'], '-dpng', '-painters');
% 
% 
close all

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
[Res,r,p]=LuPairRegressPlot(squeeze(rStim(:,1,2)),NegRate(:),ResParam);
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
print(gcf, [ResultFolder 'SepFunSingleSLMVsStimRegress.svg'], '-dsvg', '-painters');
print(gcf, [ResultFolder 'SepFunSingleSLMVsStimRegress.eps'], '-depsc', '-painters');
print(gcf, [ResultFolder 'SepFunSingleSLMVsStimRegress.tif'], '-dtiffn', '-painters');
% print(gcf, [ResultFolder 'SingleSLMVsSpeed.png'], '-dpng', '-painters');







close all
TargetOutDeg=zeros(length(TargetCellList),1);
TargetOutMinStrength=zeros(length(TargetCellList),1);
TargetOutStrength=zeros(length(TargetCellList),length(iscell));

crit_pTemp=crit_pAll;

for iCell=1:length(TargetCellList)
    if SuccTarget(iCell)
       TargetC=TargetCellList(iCell);
       temp1=statCellRes(iCell,PowerTargetI(iCell)).delta;
       temp2=statCellRes(iCell,PowerTargetI(iCell)).p>crit_pTemp;
       % temp2=[];
       temp1(temp2)=NaN;
       TargetOutDeg(iCell)=sum(~isnan(temp1));
       temp3=temp1;
       % temp3(isnan(temp3))=0;
       TargetOutMinStrength(iCell,1)=min(abs(temp3));
       TargetOutStrength(iCell,:)=temp3;

       % TargetResAmp(iCell)=TargetCellResP(iCell,PowerTargetI(iCell))
    else

    end
end













TestStepFrame=3;
ExampleCell=find(SuccAmp(:)>0.1);
ExampleCell=find(TargetOutDeg<=40&abs(TargetOutMinStrength)>0.05&SuccAmp(:)>0.2)


        Param3DMat.EdgeParam=[0.06 0.1 0.06 0.06 0.06 0.06];
        Param3DMat.CellCenterWith=1;
        Param3DMat.CellBoundaryWidth=0.5;
        Param3DMat.PlotCenter=1;
        Param3DMat.Alpha=1;
        Param3DMat.HighLightWidth=2;
Param3DMat.HighLightColor=[0 1 0];
MapPSTHparam=PSTHparam;
MapPSTHparam.PostSLMCal=TestStepFrame;


ParamTrace.PlotType=3;
ParamTrace.statisP=0;
ParamTrace.LegendShow=0;
ParamTrace.Legend=[];

NonTargetColor=[0.9 0.6 0.0];

%%ExampleCell=14;
rawPixelMapResponse=[-600 600];


close all
for jCell=1:length(ExampleCell)
    iCell=ExampleCell(jCell);
    PostSynCell=find(abs(TargetOutStrength(iCell,:))>0);
    PostSynCell=setdiff(PostSynCell,TargetCellList(iCell));

    PointTemp=AlignedInfoTable.Point(AlignedInfoTable.PointTargetCell==TargetCellList(iCell));
    PointTemp=unique(PointTemp);


    PostSLMShow=3;
    PlotCellList=[TargetCellList(iCell) PostSynCell];

    tempColor=[GroupColor(TargetCellListFunGroup(iCell),:);repmat(NonTargetColor,length(PostSynCell),1)];
    
    figure;
    t1=subplotLU(1,2,1,1,[0.1 0.1 0.1 0.05 0.04 0.04]);
    % MultiPlanes3DShow(permute(CaData.PlaneMeanImg,[2 1 3]), cellBoundary(PlotCellList), NeuronPos3D(PlotCellList,:), [], Zdepth, tempColor, [0 400],Param3DMat)
    MultiPlanes3DShow_HightLight(permute(CaData.PlaneMeanImg,[2 1 3]), [], [],[], Zdepth, tempColor, [0 400],SLMPos3D(PointTemp,:),Param3DMat)
 
    set(gca,'xlim',[0 512],'ylim',[0 512]);axis off
    axis off
    bar1=colorbar(t1);
    bar1.Location='northoutside';
    bar1.Position=[0.2 0.92 0.1 0.01];
    bar1.Label.String='F.';
    colormap(t1,bone);

    papersizePX=[0 0 16 6];

   I1=find(SLMInfoTable.PointTargetCell==TargetCellList(iCell)&abs(SLMInfoTable.UncagingLaserPower-PVpower(2))<0.1)

   I1=find(AlignedInfoTable.PointTargetCell==TargetCellList(iCell)&abs(AlignedInfoTable.UncagingLaserPower-PVpower(PowerTargetI(iCell)))<0.1)
   PSTHall = CalPSTH_suite2pBin(confSet, MapPSTHparam,AlignedInfoTable.suite2pInd(I1));
   PSTHtestCell=squeeze(nanmean(PSTHall,3));
   PSTHtestCell=SmoothDecDim3(PSTHtestCell,1);

   clear tempLabel
   tempLabel{1}='SLM';
   for il=1:length(PostSynCell)
       tempLabel{il+1}=num2str(PostSynCell(il));
   end

   tempColor=[repmat(NonTargetColor,length(PostSynCell)+1,1)];

   Param3DMat.PlotCenter=0;
   Param3DMat.HighLightColor=[0 1 0];
   Param3DMat.HighLightWidth=2;
   Param3DMat.HighLightWidth=2;

    t1=subplotLU(1,2,1,2,[0.1 0.1 0.1 0.05 0.04 0.04]);
   MultiPlanes3DShow(permute(PSTHtestCell,[1 2 3]), cellBoundary(PlotCellList), NeuronPos3D(PlotCellList,:), tempLabel, Zdepth, tempColor, [-400 400],Param3DMat)
   MultiPlanes3DShow_HightLight(permute(PSTHtestCell,[1 2 3]), cellBoundary(PostSynCell), NeuronPos3D(PostSynCell,:),[], Zdepth, tempColor, [-400 400],SLMPos3D(PointTemp,:),Param3DMat)
   set(gca,'clim',rawPixelMapResponse)
   colormap(t1,ParamNet.ResponseMap);

   axis off
   bar2=colorbar;
    bar2.Location='northoutside';
    bar2.Position=[0.6 0.92 0.1 0.01];
    bar2.Label.String='Responses';
    papersizePX=[0 0 30 15];
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
    % print(gcf, [EGTestFolder 'MapCellTarget' num2str(iCell) '.svg'], '-dsvg', '-painters');
    % print(gcf, [EGTestFolder 'MapCellTarget' num2str(iCell) '.eps'], '-depsc', '-painters');
    % saveas(gcf,[EGTestFolder 'MapCellTarget' num2str(iCell) '.eps'], 'epsc')
    saveas(gcf,[EGTestFolder 'MapCellTarget' num2str(iCell) ], 'tif')
    % 
    % print(gcf, [EGTestFolder 'MapCellTarget' num2str(iCell) '.tif'], '-dtiff', '-painters');






   [~,planeTemp]=ismember(SLMPos3D(PointTemp,3),Zdepth);
   figure;
   for iPowerLevel=1:length(PVpower)

       subplotLU(1,length(PVpower),1,iPowerLevel,[0.01 0.15 0.08 0.05 0.08 0.01])

       I1=find(AlignedInfoTable.PointTargetCell==TargetCellList(iCell)&abs(AlignedInfoTable.UncagingLaserPower-PVpower(iPowerLevel))<0.1)
       if length(I1)>1
       PSTHtemp = CalPSTH_suite2pBin(confSet, MapPSTHparam,AlignedInfoTable.suite2pInd(I1));
       PSTHtemp=squeeze(mean(PSTHtemp(:,:,:,planeTemp),3));
       roiHeat = Center2neighbour(PSTHtemp, SLMPos3D(PointTemp,1:2), SLMPosInfo.ROIparam.NeighbourHfWidthPixel);
       imagesc(SmoothDec(roiHeat',0.5));
       % hold on;contour(SmoothDec(roiHeat',0.5),[400 800],'g-','LineWidth',2)
             [SLMResTemp, mergedRegion, contourPixels] = ROIcontour_center(roiHeat,SLMPosInfo.ROIparam);
      if SLMResTemp
          hold on;
         plot(contourPixels(:,1),contourPixels(:,2),'g.');

      end

      else
       imagesc(zeros(size(roiHeat')));
      


       end
       title(['Power ' num2str(iPowerLevel)])

       colormap(ParamNet.ResponseMap);
       set(gca,'clim',[-1200 1200]);
 
       axis off

   end
   LuFontStandard 
   bar1=colorbar;
    bar1.Location='eastoutside';
    bar1.Position=[0.9 0.3 0.01 0.4];
    bar1.Label.String='Responses';
    papersizePX=[0 0 16 6];
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
    % print(gcf, [EGTestFolder 'CellTarget' num2str(iCell) '.svg'], '-dsvg', '-painters');
    print(gcf, [EGTestFolder 'ROICellTarget' num2str(iCell) '.eps'], '-depsc', '-painters');
    print(gcf, [EGTestFolder 'ROICellTarget' num2str(iCell) '.tif'], '-dtiff', '-painters');

   
     I1=find(AlignedInfoTable.PointTargetCell==TargetCellList(iCell)&abs(AlignedInfoTable.UncagingLaserPower-PVpower(PowerTargetI(iCell)))<0.1);

     tempPlotTrace1{1}=squeeze(AlignedtempNData(TargetCellList(iCell),:,I1))';

     ParamTrace.PlotType=3
     figure;
     subplotLU(1,2,1,1,[0.1 0.1 0.1 0.1 0.05 0.1])

     RateHist_GroupPlot(TimBinFrame+0.5,tempPlotTrace1,Param3DMat.HighLightColor,ParamTrace)
            hold on;
     % text(-10,0.1,['n = ' num2str(CellSampleN(iCell,iPower))])
     set(gca,'xlim',[-PSTHparam.PreSLMCal PSTHparam.PostSLMCal],'xtick',[-PSTHparam.PreSLMCal:5:PSTHparam.PostSLMCal],'ylim',[-0.3 0.4],'Box','off')
     PostSynCellAveTrace=squeeze(nanmean(AlignedtempNData(PostSynCell,:,I1),3));

     for iiPCell=1:length(PostSynCell)  
         tempPost{iiPCell}=squeeze(AlignedtempNData(PostSynCell(iiPCell),:,I1))';
     end
     plot(TimBinFrame+0.5,PostSynCellAveTrace,'Color',NonTargetColor,'LineWidth',2);
     plot(TimBinFrame+0.5,zeros(size(TimBinFrame)),'k:');
     ylabel('Normalized F.');
     xlabel('Frames from SLM')


      LuLegend([-8 -8;0.3 0.25;-7 -7;0.3 0.25],0,{'Target cell',[num2str(length(PostSynCell) ) ' impacted cells']},[Param3DMat.HighLightColor;NonTargetColor],8);


      TempAdj=zeros(length(iscell))+NaN;
      TempAdj(TargetCellList(iCell),PostSynCell)=statCellRes(iCell,PowerTargetI(iCell)).delta(PostSynCell);
      ImpactNode=nansum(TempAdj);

      TempNode=zeros(length(iscell),1);
      TempNode([TargetCellList(iCell) PostSynCell])=1;
      NodeTarget=zeros(length(iscell),1);




      t2=subplotLU(1,2,1,2,[0.1 0.1 0.1 0.1 0.05 0.1])
      [GgraphOut.G,GgraphOut.p]=MarkovState_HeatStrPlot(TempAdj,TempNode,TempAdj,ImpactNode,ParamNet.NodeColor,ParamNet.ResponseMap,ParamNet.ResponseLim);
      GgraphOut.p.NodeColor=ParamNet.NodeColor;
      GgraphOut.p.NodeColor(TargetCellList(iCell),:)=Param3DMat.HighLightColor;
      GgraphOut.p.NodeColor(PostSynCell,:)=repmat(NonTargetColor,length(PostSynCell),1);


radius=5;
theta=ClockWiseGraph(GgraphOut.G,GgraphOut.p,radius);
% p.MarkerSize=p.MarkerSize;
radius
CircosBand(1).rband=[radius+0.5;radius+2];
CircosBand(1).Values=ImpactNode;
CircosBand(1).Clim=ParamNet.ResponseLim;
CircosBand(1).theta=theta;
CircosBand(1).Colormap=ParamNet.ResponseMap;
CircosBand(1).arc_gap=abs(theta(2)-theta(1))/10;
AddCircosBand(CircosBand);
GgraphOut.p.ArrowSize=10;

colormap(ParamNet.ResponseMap);
set(gca,'clim',ParamNet.ResponseLim);
bar2=colorbar(t2);
bar2.Location='eastoutside';
bar2.Position=[0.92 0.2 0.02 0.6];
bar2.Label.String='Response';
bar2.Ticks=union(ParamNet.ResponseLim,0);

    papersizePX=[0 0 17.5 8];
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
    print(gcf, [EGTestFolder 'FC' num2str(iCell) '.svg'], '-dsvg', '-painters');
    % print(gcf, [EGTestFolder 'FC' num2str(iCell) '.eps'], '-depsc', '-painters');
    print(gcf, [EGTestFolder 'FC' num2str(iCell) '.tif'], '-dtiff', '-painters');

    






     % RateHist_GroupPlot(TimBinFrame+0.5,tempPost,tempColor,ParamTrace)

      %   CellSampleN(iCell,iPower)=length(I1);
      % 
      %   if length(I1)>1
      %      % CellResponse{iCell,iPower}=squeeze(AlignedtempNData(:,:,I1));
      %      % TargetResponse{iCell,iPower}=squeeze(AlignedtempNData(TargetCellList(iCell),:,I1))';
      % 
      %   % elseif length(I1)>1
      %      CellResponse{iCell,iPower}=nanmean(AlignedtempNData(:,:,I1),3);
      %      TargetResponse{iCell,iPower}=squeeze(AlignedtempNData(TargetCellList(iCell),:,I1))';
      %   % else
      %      preSLMdata=squeeze(AlignedtempNData(:,1:PSTHparam.PreSLMCal,I1));
      %      postSLMdata=squeeze(AlignedtempNData(:,PSTHparam.PreSLMCal+[1:PostSLMShow],I1));
      %      for jCell=1:length(iscell)
      %          temp1=preSLMdata(jCell,:,:);
      %          temp2=postSLMdata(jCell,:,:);
      %      end
      %      statCellRes(iCell,iPower).p=p;
      %      statCellRes(iCell,iPower).t=t;
      %      statCellRes(iCell,iPower).delta=temp3;
      % 
      %      TargetCellResP(iCell,iPower)=p(TargetCellList(iCell));
      %      TargetCellResR(iCell,iPower)=temp3(TargetCellList(iCell));
      %      clear p t;
      %   end
      % 
close all

      % 

end



close all
for jCell=1:length(ExampleCell)
    iCell=ExampleCell(jCell);
    PostSynCell=find(abs(TargetOutStrength(iCell,:))>0);
    PostSynCell=setdiff(PostSynCell,TargetCellList(iCell));

    PointTemp=AlignedInfoTable.Point(AlignedInfoTable.PointTargetCell==TargetCellList(iCell));
    PointTemp=unique(PointTemp);


    PostSLMShow=3;
    PlotCellList=[TargetCellList(iCell) PostSynCell];

    tempColor=[GroupColor(TargetCellListFunGroup(iCell),:);repmat(NonTargetColor,length(PostSynCell),1)];
    

   I1=find(SLMInfoTable.PointTargetCell==TargetCellList(iCell)&abs(SLMInfoTable.UncagingLaserPower-PVpower(2))<0.1)

   I1=find(AlignedInfoTable.PointTargetCell==TargetCellList(iCell)&abs(AlignedInfoTable.UncagingLaserPower-PVpower(PowerTargetI(iCell)))<0.1)
   PSTHall = CalPSTH_suite2pBin(confSet, MapPSTHparam,AlignedInfoTable.suite2pInd(I1));
   PSTHtestCell=squeeze(nanmean(PSTHall,3));
   PSTHtestCell=SmoothDecDim3(PSTHtestCell,1);



   [~,planeTemp]=ismember(SLMPos3D(PointTemp,3),Zdepth);
   figure;
   for iPowerLevel=1:length(PVpower)

       subplotLU(length(PVpower),1,iPowerLevel,1,[0.01 0.15 0.08 0.05 0.08 0.04])

       I1=find(AlignedInfoTable.PointTargetCell==TargetCellList(iCell)&abs(AlignedInfoTable.UncagingLaserPower-PVpower(iPowerLevel))<0.1)
       if length(I1)>1
       PSTHtemp = CalPSTH_suite2pBin(confSet, MapPSTHparam,AlignedInfoTable.suite2pInd(I1));
       PSTHtemp=squeeze(mean(PSTHtemp(:,:,:,planeTemp),3));
       roiHeat = Center2neighbour(PSTHtemp, SLMPos3D(PointTemp,1:2), SLMPosInfo.ROIparam.NeighbourHfWidthPixel);
       imagesc(SmoothDec(roiHeat',0.5));
       % hold on;contour(SmoothDec(roiHeat',0.5),[400 800],'g-','LineWidth',2)
             [SLMResTemp, mergedRegion, contourPixels] = ROIcontour_center(roiHeat,SLMPosInfo.ROIparam);
      % % if SLMResTemp
      % %     hold on;
      % %    plot(contourPixels(:,1),contourPixels(:,2),'g.');
      % % 
      % % end

      else
       imagesc(zeros(size(roiHeat')));
      


       end
       % title(['Power ' num2str(iPowerLevel)])

       colormap(ParamNet.ResponseMap);
       set(gca,'clim',rawPixelMapResponse);
 
       axis off

   end
   LuFontStandard 
   % bar1=colorbar;
   %  bar1.Location='eastoutside';
   %  bar1.Position=[0.9 0.3 0.01 0.4];
   %  bar1.Label.String='Responses';
    papersizePX=[0 0 5 17];
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
    print(gcf, [EGTestFolder 'ROICellTarget' num2str(iCell) 'Vertical.eps'], '-depsc', '-painters');
    print(gcf, [EGTestFolder 'ROICellTarget' num2str(iCell) 'Vertical.tif'], '-dtiff', '-painters');






     % RateHist_GroupPlot(TimBinFrame+0.5,tempPost,tempColor,ParamTrace)

      %   CellSampleN(iCell,iPower)=length(I1);
      % 
      %   if length(I1)>1
      %      % CellResponse{iCell,iPower}=squeeze(AlignedtempNData(:,:,I1));
      %      % TargetResponse{iCell,iPower}=squeeze(AlignedtempNData(TargetCellList(iCell),:,I1))';
      % 
      %   % elseif length(I1)>1
      %      CellResponse{iCell,iPower}=nanmean(AlignedtempNData(:,:,I1),3);
      %      TargetResponse{iCell,iPower}=squeeze(AlignedtempNData(TargetCellList(iCell),:,I1))';
      %   % else
      %      preSLMdata=squeeze(AlignedtempNData(:,1:PSTHparam.PreSLMCal,I1));
      %      postSLMdata=squeeze(AlignedtempNData(:,PSTHparam.PreSLMCal+[1:PostSLMShow],I1));
      %      for jCell=1:length(iscell)
      %          temp1=preSLMdata(jCell,:,:);
      %          temp2=postSLMdata(jCell,:,:);
      %      end
      %      statCellRes(iCell,iPower).p=p;
      %      statCellRes(iCell,iPower).t=t;
      %      statCellRes(iCell,iPower).delta=temp3;
      % 
      %      TargetCellResP(iCell,iPower)=p(TargetCellList(iCell));
      %      TargetCellResR(iCell,iPower)=temp3(TargetCellList(iCell));
      %      clear p t;
      %   end
      % 


      % 

end


rawPixelMapResponse=[-400 400];
close all
pEdge=[0.1 0.05 0.15 0.1 0.05 0.1];
PosMatX=[0.3 0.3 0.35];
PosMatY=[0.75 0.75 0.75];
pEdge=[0.1 0.05 0.15 0.1];


for jCell=1:length(ExampleCell)
    iCell=ExampleCell(jCell);
    PostSynCell=find(abs(TargetOutStrength(iCell,:))>0);
    PostSynCell=setdiff(PostSynCell,TargetCellList(iCell));

    PointTemp=AlignedInfoTable.Point(AlignedInfoTable.PointTargetCell==TargetCellList(iCell));
    PointTemp=unique(PointTemp);


    PostSLMShow=3;
    PlotCellList=[TargetCellList(iCell) PostSynCell];

    tempColor=[GroupColor(TargetCellListFunGroup(iCell),:);repmat(NonTargetColor,length(PostSynCell),1)];
    

    I1=find(AlignedInfoTable.PointTargetCell==TargetCellList(iCell)&abs(AlignedInfoTable.UncagingLaserPower-PVpower(PowerTargetI(iCell)))<0.1);

    tempPlotTrace1{1}=squeeze(AlignedtempNData(TargetCellList(iCell),:,I1))';

    ParamTrace.PlotType=3;
    figure;
    subplotPosLu(PosMatX,PosMatY,1,1,pEdge);

    RateHist_GroupPlot(TimBinFrame+0.5,tempPlotTrace1,Param3DMat.HighLightColor,ParamTrace)
    hold on;
     % text(-10,0.1,['n = ' num2str(CellSampleN(iCell,iPower))])
     set(gca,'xlim',[-PSTHparam.PreSLMCal PSTHparam.PostSLMCal],'xtick',[-PSTHparam.PreSLMCal:5:PSTHparam.PostSLMCal],'ylim',[-0.3 0.4],'Box','off')
     PostSynCellAveTrace=squeeze(nanmean(AlignedtempNData(PostSynCell,:,I1),3));

     for iiPCell=1:length(PostSynCell)  
         tempPost{iiPCell}=squeeze(AlignedtempNData(PostSynCell(iiPCell),:,I1))';
     end
     plot(TimBinFrame+0.5,PostSynCellAveTrace,'Color',NonTargetColor,'LineWidth',2);
     plot(TimBinFrame+0.5,zeros(size(TimBinFrame)),'k:');
     ylabel('Normalized F.');
     xlabel('Frames from SLM')


     LuLegend([-8 -8;0.3 0.25;-7 -7;0.3 0.25],0,{'Target cell',[num2str(length(PostSynCell) ) ' impacted cells']},[Param3DMat.HighLightColor;NonTargetColor],8);



   % t2=subplotLU(1,3,1,2,pEdge)
    t2=subplotPosLu(PosMatX,PosMatY,1,2,pEdge);

   I1=find(AlignedInfoTable.PointTargetCell==TargetCellList(iCell)&abs(AlignedInfoTable.UncagingLaserPower-PVpower(PowerTargetI(iCell)))<0.1)
   PSTHall = CalPSTH_suite2pBin(confSet, MapPSTHparam,AlignedInfoTable.suite2pInd(I1));
   PSTHtestCell=squeeze(nanmean(PSTHall,3));
   PSTHtestCell=SmoothDecDim3(PSTHtestCell,1);

   Param3DMatAlpha=Param3DMat;
   Param3DMatAlpha.Alpha=1;
   % MultiPlanes3DShow(permute(PSTHtestCell,[1 2 3]), [], [], [], Zdepth, tempColor, [-400 400],Param3DMatAlpha)
   MultiPlanes3DShow_HightLight(permute(PSTHtestCell,[1 2 3]), [], [],[], Zdepth, tempColor, rawPixelMapResponse,SLMPos3D(PointTemp,:),Param3DMatAlpha)
   set(gca,'clim',rawPixelMapResponse)
   colormap(t2,ParamNet.ResponseMap);
% bar3=colorbar(t3);
% bar3.Location='northoutside';
% bar3.Position=[0.8 0.92 0.1 0.03];
% bar3.Label.String='Dist. from Target';
% bar3.Ticks=union(ParamNetDist.ScoreLim,0);
% colormap(ParamNet.ResponseMap);
set(gca,'clim',rawPixelMapResponse);

axis off
      TempAdjSub=zeros(length(PostSynCell)+1)+NaN;
      TempAdjSub(1,2:end)=statCellRes(iCell,PowerTargetI(iCell)).delta(PostSynCell);
      ImpactNodeSub=nansum(TempAdjSub);

      TempNodeSub=zeros(length(PostSynCell)+1,1)+1;
      % TempNodeSub([TargetCellList(iCell) PostSynCell])=1;
      NodeTarget=zeros(length(iscell),1);


      TempAdj=zeros(length(iscell))+NaN;
      TempAdj(TargetCellList(iCell),PostSynCell)=statCellRes(iCell,PowerTargetI(iCell)).delta(PostSynCell);
      ImpactNode=nansum(TempAdj);

      TempNode=zeros(length(iscell),1);
      TempNode([TargetCellList(iCell) PostSynCell])=1;
      NodeTarget=zeros(length(iscell),1);




   hold on;

    axis off

% colormap(ParamNet.ResponseMap);
% set(gca,'clim',ParamNet.ResponseLim);




bar2=colorbar(t2);
bar2.Location='northoutside';
bar2.Position=[0.4 0.88 0.1 0.01];
bar2.Label.String='Response Raw';
bar2.Ticks=union(rawPixelMapResponse,0);


 %    t3=subplotLU(1,3,1,3,pEdge)
 % 
 %    plot(1,1);
 %   set(gca, 'Visible', 'off');
 % axis off


% bar3=colorbar(t3);
% bar3.Location='northoutside';
% bar3.Position=[0.8 0.92 0.1 0.03];
% bar3.Label.String='Response Norm.';
% bar3.Ticks=union(ParamNet.ResponseLim,0);
% axis off

    papersizePX=[0 0 18 9];
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
    print(gcf, [EGTestFolder 'FC' num2str(iCell) 'NoLinks.tif'], '-dtiff', '-painters');

    % print(gcf, [EGTestFolder 'FC' num2str(iCell) 'NoLinks.svg'], '-dsvg', '-painters');
    % print(gcf, [EGTestFolder 'FC' num2str(iCell) 'NoLinks.eps'], '-depsc', '-painters');
% 
% 
% 
% close all

    % % bar2=colorbar;
    % % bar2.Location='northoutside';
    % % bar2.Position=[0.6 0.92 0.1 0.01];
    % % bar2.Label.String='Responses';
    % % papersizePX=[0 0 30 15];


end


close all
for jCell=1:length(ExampleCell)
    iCell=ExampleCell(jCell);
    PostSynCell=find(abs(TargetOutStrength(iCell,:))>0);
    PostSynCell=setdiff(PostSynCell,TargetCellList(iCell));

    PointTemp=AlignedInfoTable.Point(AlignedInfoTable.PointTargetCell==TargetCellList(iCell));
    PointTemp=unique(PointTemp);


    PostSLMShow=3;
    PlotCellList=[TargetCellList(iCell) PostSynCell];

    tempColor=[GroupColor(TargetCellListFunGroup(iCell),:);repmat(NonTargetColor,length(PostSynCell),1)];
    

    I1=find(AlignedInfoTable.PointTargetCell==TargetCellList(iCell)&abs(AlignedInfoTable.UncagingLaserPower-PVpower(PowerTargetI(iCell)))<0.1);

    tempPlotTrace1{1}=squeeze(AlignedtempNData(TargetCellList(iCell),:,I1))';

    ParamTrace.PlotType=3;
    figure;
    subplotPosLu(PosMatX,PosMatY,1,1,pEdge);

    RateHist_GroupPlot(TimBinFrame+0.5,tempPlotTrace1,Param3DMat.HighLightColor,ParamTrace)
    hold on;
     % text(-10,0.1,['n = ' num2str(CellSampleN(iCell,iPower))])
     set(gca,'xlim',[-PSTHparam.PreSLMCal PSTHparam.PostSLMCal],'xtick',[-PSTHparam.PreSLMCal:5:PSTHparam.PostSLMCal],'ylim',[-0.3 0.4],'Box','off')
     PostSynCellAveTrace=squeeze(nanmean(AlignedtempNData(PostSynCell,:,I1),3));

     for iiPCell=1:length(PostSynCell)  
         tempPost{iiPCell}=squeeze(AlignedtempNData(PostSynCell(iiPCell),:,I1))';
     end
     plot(TimBinFrame+0.5,PostSynCellAveTrace,'Color',NonTargetColor,'LineWidth',2);
     plot(TimBinFrame+0.5,zeros(size(TimBinFrame)),'k:');
     ylabel('Normalized F.');
     xlabel('Frames from SLM')


      LuLegend([-8 -8;0.3 0.25;-7 -7;0.3 0.25],0,{'Target cell',[num2str(length(PostSynCell) ) ' impacted cells']},[Param3DMat.HighLightColor;NonTargetColor],8);



   % t2=subplotLU(1,3,1,2,pEdge)
   t2=subplotPosLu(PosMatX,PosMatY,1,2,pEdge);

   I1=find(AlignedInfoTable.PointTargetCell==TargetCellList(iCell)&abs(AlignedInfoTable.UncagingLaserPower-PVpower(PowerTargetI(iCell)))<0.1)
   PSTHall = CalPSTH_suite2pBin(confSet, MapPSTHparam,AlignedInfoTable.suite2pInd(I1));
   PSTHtestCell=squeeze(nanmean(PSTHall,3));
   PSTHtestCell=SmoothDecDim3(PSTHtestCell,1);

   Param3DMatAlpha=Param3DMat;
   Param3DMatAlpha.Alpha=0.8
   % MultiPlanes3DShow(permute(PSTHtestCell,[1 2 3]), [], [], [], Zdepth, tempColor, [-400 400],Param3DMatAlpha)
   MultiPlanes3DShow_HightLight(permute(PSTHtestCell,[1 2 3]), [], [],[], Zdepth, tempColor, rawPixelMapResponse,SLMPos3D(PointTemp,:),Param3DMatAlpha)
   set(gca,'clim',rawPixelMapResponse)
   colormap(t2,ParamNet.ResponseMap);
% bar3=colorbar(t3);
% bar3.Location='northoutside';
% bar3.Position=[0.8 0.92 0.1 0.03];
% bar3.Label.String='Dist. from Target';
% bar3.Ticks=union(ParamNetDist.ScoreLim,0);
% colormap(ParamNet.ResponseMap);
set(gca,'clim',rawPixelMapResponse);

axis off
      axis off
      TempAdjSub=zeros(length(PostSynCell)+1)+NaN;
      TempAdjSub(1,2:end)=statCellRes(iCell,PowerTargetI(iCell)).delta(PostSynCell);
      ImpactNodeSub=nansum(TempAdjSub);

      TempNodeSub=zeros(length(PostSynCell)+1,1)+1;
      % TempNodeSub([TargetCellList(iCell) PostSynCell])=1;
      NodeTarget=zeros(length(iscell),1);


      TempAdj=zeros(length(iscell))+NaN;
      TempAdj(TargetCellList(iCell),PostSynCell)=statCellRes(iCell,PowerTargetI(iCell)).delta(PostSynCell);
      ImpactNode=nansum(TempAdj);

      TempNode=zeros(length(iscell),1);
      TempNode([TargetCellList(iCell) PostSynCell])=1;
      NodeTarget=zeros(length(iscell),1);




   hold on;


    % theta=ClockWiseGraph(GgraphOut.G,GgraphOut.p,radius);
    [GgraphOut.G,GgraphOut.p]=MarkovState_HeatStrPlot(TempAdj,TempNode,TempAdj,ImpactNode,ParamNet.NodeColor,ParamNet.ResponseMap,ParamNet.ResponseLim);

    GgraphOut.p.XData = NeuronPos3D(:,2)';
    GgraphOut.p.YData = NeuronPos3D(:,1)';
    GgraphOut.p.ZData = NeuronPos3D(:,3)';
    GgraphOut.p.NodeLabel=[];
    GgraphOut.p.ArrowSize=8;
    GgraphOut.p.Marker='none';
% colormap(ParamNet.ResponseMap);
% set(gca,'clim',ParamNet.ResponseLim);
% bar2=colorbar(t2);
% bar2.Location='northoutside';
% bar2.Position=[0.6 0.92 0.1 0.03];
% bar2.Label.String='Response';
% bar2.Ticks=union([-800 800],0);
    axis off


      
      TempAdj=zeros(length(iscell))+NaN;
      TempAdj(TargetCellList(iCell),PostSynCell)=statCellRes(iCell,PowerTargetI(iCell)).delta(PostSynCell);
      ImpactNode=nansum(TempAdj);

      TempNode=zeros(length(iscell),1);
      TempNode([TargetCellList(iCell) PostSynCell])=1;
      NodeTarget=zeros(length(iscell),1);




      t3=subplotPosLu(PosMatX,PosMatY,1,3,pEdge);
      [GgraphOut.G,GgraphOut.p]=MarkovState_HeatStrPlot(TempAdj,TempNode,TempAdj,ImpactNode,ParamNet.NodeColor,ParamNet.ResponseMap,ParamNet.ResponseLim);
      GgraphOut.p.NodeColor=ParamNet.NodeColor;
      GgraphOut.p.NodeColor(TargetCellList(iCell),:)=Param3DMat.HighLightColor;
      GgraphOut.p.NodeColor(PostSynCell,:)=repmat(NonTargetColor,length(PostSynCell),1);


radius=5;
theta=ClockWiseGraph(GgraphOut.G,GgraphOut.p,radius);
% p.MarkerSize=p.MarkerSize;
radius
clear CircosBand
CircosBand(1).rband=[radius+0.5;radius+2];
CircosBand(1).Values=ImpactNode;
CircosBand(1).Clim=ParamNet.ResponseLim;
CircosBand(1).theta=theta;
CircosBand(1).Colormap=ParamNet.ResponseMap;
CircosBand(1).arc_gap=abs(theta(2)-theta(1))/10;

% CircosBand(end+1).rband=CircosBand(end).rband+1.2;
% CircosBand(end).Values=CellDist(TargetCellList(iCell),:);
% CircosBand(end).Clim=ParamNetDist.ScoreLim;
% CircosBand(end).theta=theta;
% CircosBand(end).Colormap=ParamNetDist.ScoreMap;
% CircosBand(end).arc_gap=abs(theta(2)-theta(1))/10;

AddCircosBand(CircosBand);
GgraphOut.p.ArrowSize=10;



colormap(t3,ParamNet.ResponseMap);
% bar3=colorbar(t3);
% bar3.Location='northoutside';
% bar3.Position=[0.8 0.92 0.1 0.03];
% bar3.Label.String='Dist. from Target';
% bar3.Ticks=union(ParamNetDist.ScoreLim,0);
% colormap(ParamNet.ResponseMap);
set(gca,'clim',ParamNet.ResponseLim);


    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));



bar2=colorbar(t2);
bar2.Location='northoutside';
bar2.Position=[0.4 0.88 0.1 0.01];
bar2.Label.String='Response Raw';
bar2.Ticks=union(rawPixelMapResponse,0);




bar3=colorbar(t3);
bar3.Location='northoutside';
bar3.Position=[0.7 0.88 0.1 0.01];

bar3.Label.String='Response Norm.';
bar3.Ticks=union(ParamNet.ResponseLim,0);
axis off


    papersizePX=[0 0 22 10];
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
    print(gcf, [EGTestFolder 'FC' num2str(iCell) 'Links.tif'], '-dtiff', '-painters');
    % 
    % % print(gcf, [EGTestFolder 'FC' num2str(iCell) 'NoLinks.svg'], '-dsvg', '-painters');
    % print(gcf, [EGTestFolder 'FC' num2str(iCell) 'Links.eps'], '-depsc', '-painters');
% 
% 
% 
% close all

    % % bar2=colorbar;
    % % bar2.Location='northoutside';
    % % bar2.Position=[0.6 0.92 0.1 0.01];
    % % bar2.Label.String='Responses';
    % % papersizePX=[0 0 30 15];


end

