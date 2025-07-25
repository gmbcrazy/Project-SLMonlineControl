function OfflineSLMAbsSpeedControl_OneFOV(FOVtemp, Suite2pDataKeywords, suite2pFOVPathLocalTemp,PSTHparam,ResultFolderSpeed)


% GroupLabel={'L','S','N'};
% nGroup=length(GroupLabel);
% GroupColor=[255 51 153;91 20 212;121 247 111]/255;
% NodeColor=repmat([0.9 0.9 0.9],length(iscell),1);
% 
papersizePX=[0 0 33 9];
papersizePXX=[0 0 22 16];
TestStepFrame=PSTHparam.TestStepFrame;
% suite2pFOVPathLocalTemp=suite2pFOVPathLocal{iFOV};
resultPaths = findAllFoldersKeyWords(suite2pFOVPathLocalTemp, Suite2pDataKeywords);



metaDataPath=Get_ExpDataFolder(resultPaths{1},'Step1Basic',{'Step1Meta.mat','BehAll.mat','.png'});
load([metaDataPath 'Step1Meta.mat'])
confSet=SLMPosInfo.confSetFinal
rSpeed=CorrResults.rSpeed;
rStim=CorrResults.rStim;

Zdepth = confSet.scan_Z + confSet.ETL;
suite2pPath = findAllFoldersKeyWords(suite2pFOVPathLocalTemp, Suite2pDataKeywords);
[~,Session,~]=fileparts(suite2pFOVPathLocalTemp(1:end-1));
SaveP1=[ResultFolderSpeed Session '\'];
mkdir(SaveP1);

confSet.save_path0=suite2pPath{1};
[NeuronPos3D,NeuronPos3DRaw,CaData,CaDataPlane,Neuronstat]=Extract_Suite2p(confSet);
% load('C:\Users\zhangl33\Projects\Project-SLMonlineControl\CodeAndPackages\subfun\Color\colorMapPN3.mat','colorMapPN1')
[~, ~, CellMedCenter, cellBoundary, ~] = Suite2pCellIDMapFromStat(CaData.statCell, [confSet.SLM_Pixels_Y confSet.SLM_Pixels_X]);
Cell3DPos=[CellMedCenter Zdepth(CaData.CellPlaneID)'];
CellDist=squareform(pdist(Cell3DPos));

iData=PSTHparam.iData;
PSTHparam.TestStepFrame=3;
crit_pAll=PSTHparam.pTh;

iscell=find(CaData.iscell(CaData.iscell(:,1)>0));
PVpower=xmlPower2PVpower(SLMTestInfo.confSet.UncagingLaserPower);
PVpower=intersect(round(PVpower),SLMInfoTable.UncagingLaserPower);
GroupLabel={'L','S','N'};
GroupList=[1 2 3];

nGroup=length(GroupLabel);
GroupColor=[255 51 153;91 20 212;121 247 111]/255;
NodeColor=repmat([0.9 0.9 0.9],length(iscell),1);
PowerZeroColor=[0.5 0.5 0.5];
PowerZero=[0 1];
PowerZeroLabel={'SLM','FakeSLM'};
VolOut=[0 1];
VolOutLabel={'NoWhisk','Whisk'}
AwakeState=[1];
AwakeStateLabel={'Awake'}
GroupMetaName=[GroupLabel {'FakeSLM'}];
GroupMetaColor=[GroupColor;PowerZeroColor];

GroupTargetCellAll=[];
for iFun=1:length(GroupList)
    GroupTargetCellAll=[GroupTargetCellAll;[GroupTargetCell{iFun}(:) zeros(size(GroupTargetCell{iFun}(:)))+iFun]];
end
GroupTargetCellMeta=[GroupTargetCell {[]}];


% [AlignedtempNData,AlignedInfoTable,statCellRes,TargetResponse,TargetCellResP,TargetCellResR,CellSampleN]=Aligned_FromSuite2p(NData{iData},TargetCellList,iscell,Suite2pTable,PVpower,PSTHparam);
[AlignedtempNData,AlignedInfoTable,CellResponse,statCellRes,TargetResponse,TargetCellResP,TargetCellResR,CellSampleN]=Aligned_FromSuite2p(NData{iData},TargetCellList,Suite2pTable,PVpower,PSTHparam);

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

figure;

   GroupPair.CorrName='fdr';
   GroupPair.Q=0.1;
   GroupPair.Pair=[1 3 5 7;2 4 6 8];
   GroupPair.SignY=5;
   GroupPair.Plot=1;
   GroupPair.Std=1;      %%%%%%%%%using standard deviation as errorbar
   GroupPair.SamplePlot=1; %%%%%%%%%Plot Individual Sample Point
   GroupPair.SamplePairedPlot=2; %%%%%%%%%Dash line for paired comparison sample
   GroupPair.LimY=[0 6];
   GroupfromMetaIndex=[1 1 2 2 3 3 4 4];
   % GroupPair.GroupName=GroupMetaName;
   clear SpeedChangeStats

for iAwake=1:length(AwakeState)
    % for iPower=1:length(PowerZero)

        for iVol=1:length(VolOut)
            clear GroupResponse GroupSpeed statGroupRes GroupSampleN

            if iAwake==1
               subplotLU(1,length(VolOutLabel),1,iVol);
               hold on;
               xData={};
                for iFun = 1:length(GroupList)
                        I1=find(AlignedInfoTable.Group==GroupList(iFun)&AlignedInfoTable.PowerZero==0&AlignedInfoTable.VolOut==VolOut(iVol)&AlignedInfoTable.AwakeState==AwakeState(iAwake));
                        % GroupSampleN(iFun)=length(I1);
                           % preSLMSpeed=mean(squeeze(AlignedSpeed(PSTHparam.PreSLMCal+[-2:1:0],I1)),1);

                           preSLMSpeed=mean(squeeze(AlignedSpeed(1:PSTHparam.PreSLMCal,I1)),1);
                           postSLMSpeed=mean(squeeze(AlignedSpeed(PSTHparam.PreSLMCal+[1:TestStepFrame],I1)),1);
                           xData{end+1}=preSLMSpeed;
                           xData{end+1}=postSLMSpeed;
                end





                %%Add Zero power as addtional group
                        I1=find(AlignedInfoTable.PowerZero==1&AlignedInfoTable.VolOut==VolOut(iVol)&AlignedInfoTable.AwakeState==AwakeState(iAwake));
                        % preSLMSpeed=mean(squeeze(AlignedSpeed(PSTHparam.PreSLMCal+[-2:1:0],I1)),1);

                        preSLMSpeed=mean(squeeze(AlignedSpeed(1:PSTHparam.PreSLMCal,I1)),1);
                        postSLMSpeed=mean(squeeze(AlignedSpeed(PSTHparam.PreSLMCal+[1:TestStepFrame],I1)),1);

                       xData{end+1}=preSLMSpeed;
                       xData{end+1}=postSLMSpeed;
                        SpeedDecRatio(iAwake,iVol,iFun+1)=length(I1)/length(preSLMSpeed);
                        % scatter(preSLMSpeed,postSLMSpeed,'MarkerFaceColor',GroupMetaColor(iFun+1,:));
                        % 
                        % plot([0.001 10],[0.001 10],'k:');
                        % set(gca,'xlim',[0 5],'ylim',[0 5],'xscale','log','yscale','log')
                 % SpeedChangeStats=ErrorBoxPlotLU(1:length(xData),xData,GroupMetaColor([ 1 1 2 2 3 3 4 4],:),[],GroupPair,[1 1 2 2 3 3 4 4]);

                 SpeedChangeStats{iVol}=ErrorBarPlotLU(1:length(xData),xData,[],GroupMetaColor(GroupfromMetaIndex,:),2,1,[SaveP1 VolOutLabel{iVol} '.txt'],GroupPair,GroupfromMetaIndex);

                 set(gca,'ylim',GroupPair.LimY,'xtick',[2:2:8]-0.5,'XTickLabel',GroupMetaName);
                 if iVol==1
                 ylabel('Speed');
                 end


            end
        end
    % end
end

                        papersizePX=[0 0 20 8];
                        set(gcf, 'PaperUnits', 'centimeters');
                       set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));

                        % print(gcf, [SaveP1 'SpeedChange.svg'], '-dsvg', '-painters');
                        % print(gcf, [SaveP1 'SpeedChange.eps'], '-depsc', '-painters');
                        print(gcf, [SaveP1 'SpeedChange.tif'], '-dtiffn', '-painters');
                
                            close all



for iSpeedTh=1:length(PostPreDiffSpeedTh)
    ResultFolder=[SaveP1 num2str(PostPreDiffSpeedTh(iSpeedTh)) '\'];
    mkdir(ResultFolder);

    for iAwake=1:length(AwakeState)
        % for iPower=1:length(PowerZero)
            for iVol=1:length(VolOut)
                TempResultFolder=[ResultFolder AwakeStateLabel{iAwake} VolOutLabel{iVol} '\'];
                mkdir(TempResultFolder);
                clear GroupResponse GroupSpeed statGroupRes GroupSampleN
    
                    for iFun = 1:length(GroupList)
                            I1=find(AlignedInfoTable.Group==GroupList(iFun)&AlignedInfoTable.PowerZero==0&AlignedInfoTable.VolOut==VolOut(iVol)&AlignedInfoTable.AwakeState==AwakeState(iAwake));
                            % GroupSampleN(iFun)=length(I1);
    
                               preSLMSpeed=mean(squeeze(AlignedSpeed(1:PSTHparam.PreSLMCal,I1)),1);
                               postSLMSpeed=mean(squeeze(AlignedSpeed(PSTHparam.PreSLMCal+[1:TestStepFrame],I1)),1);
    
                               if iAwake==1
                               I1=I1(find(abs(postSLMSpeed-preSLMSpeed)<=PostPreDiffSpeedTh(iSpeedTh)));
                               end
                               % temp=find((postSLMSpeed-preSLMSpeed)>=0))
                               % if length(temp)<=1
                               % else
                               %    I1=
                               % end
    
                               GroupSampleN(iFun)=length(I1);
    
                            if length(I1)==1
                               GroupResponse{iFun}=squeeze(AlignedtempNData(:,:,I1));
                               GroupSpeed{iFun}=AlignedSpeed(:,I1)';
                               preSLMdata=squeeze(AlignedtempNData(:,1:PSTHparam.PreSLMCal,I1));
                               postSLMdata=squeeze(AlignedtempNData(:,PSTHparam.PreSLMCal+[1:TestStepFrame],I1));
                               for jCell=1:length(iscell)
                                   temp1=preSLMdata(jCell,:,:);
                                   temp2=postSLMdata(jCell,:,:);
                                   [~,p(jCell,1),~,t(jCell)]=ttest2(temp2(:),temp1(:));
                                   temp3(jCell)=mean(temp2(:))-mean(temp1(:));
                               end
                            elseif length(I1)>1
                               GroupResponse{iFun}=nanmean(AlignedtempNData(:,:,I1),3);           
                               GroupSpeed{iFun}=AlignedSpeed(:,I1)';
                               preSLMdata=squeeze(AlignedtempNData(:,1:PSTHparam.PreSLMCal,I1));
                               postSLMdata=squeeze(AlignedtempNData(:,PSTHparam.PreSLMCal+[1:TestStepFrame],I1));
                               for jCell=1:length(iscell)
                                   temp1=preSLMdata(jCell,:,:);
                                   temp2=postSLMdata(jCell,:,:);
                                   [~,p(jCell,1),~,t(jCell)]=ttest2(temp2(:),temp1(:));
                                   temp3(jCell)=mean(temp2(:))-mean(temp1(:));
                               end
                               statGroupRes(iFun).p=p;
                               statGroupRes(iFun).t=t;
                               statGroupRes(iFun).delta=temp3;
                    
                               clear p t temp3;
                    
                    
                            else
                               statGroupRes(iFun).p=zeros(size(iscell))+1;
                               statGroupRes(iFun).delta=zeros(size(iscell));
    
                            end
                    
                    end
    
    
    
                
    
                    %%Add Zero power as addtional group
                            I1=find(AlignedInfoTable.PowerZero==1&AlignedInfoTable.VolOut==VolOut(iVol)&AlignedInfoTable.AwakeState==AwakeState(iAwake));
                            preSLMSpeed=mean(squeeze(AlignedSpeed(1:PSTHparam.PreSLMCal,I1)),1);
                            postSLMSpeed=mean(squeeze(AlignedSpeed(PSTHparam.PreSLMCal+[1:TestStepFrame],I1)),1);
                            if iAwake==1
                               I1=I1(find(abs(postSLMSpeed-preSLMSpeed)<=PostPreDiffSpeedTh(iSpeedTh)));
                            end
    
                            GroupSampleN(iFun+1)=length(I1);
                            if length(I1)==1
                               GroupResponse{iFun+1}=squeeze(AlignedtempNData(:,:,I1));
                               GroupSpeed{iFun+1}=AlignedSpeed(:,I1)';
                            elseif length(I1)>1
                               GroupResponse{iFun+1}=nanmean(AlignedtempNData(:,:,I1),3);           
                               GroupSpeed{iFun+1}=AlignedSpeed(:,I1)';
                               preSLMdata=squeeze(AlignedtempNData(:,1:PSTHparam.PreSLMCal,I1));
                               postSLMdata=squeeze(AlignedtempNData(:,PSTHparam.PreSLMCal+[1:TestStepFrame],I1));
                               for jCell=1:length(iscell)
                                   temp1=preSLMdata(jCell,:,:);
                                   temp2=postSLMdata(jCell,:,:);
                                   [~,p(jCell,1),~,t(jCell)]=ttest2(temp2(:),temp1(:));
                                   temp3(jCell)=mean(temp2(:))-mean(temp1(:));
                               end
                               statGroupRes(iFun+1).p=p;
                               statGroupRes(iFun+1).t=t;
                               statGroupRes(iFun+1).delta=temp3;
                    
                               clear p t temp3;
                    
                    
                            else
                               statGroupRes(iFun+1).p=zeros(size(iscell))+1;
                               statGroupRes(iFun+1).delta=zeros(size(iscell));
    
                            end
                     %%Add Zero power as addtional group
    
                        figure;
                        % TargetCellList(iCell)
                        % subplotLU(length(TargetCellList),length(PVpower),iCell,iPower,P)
                        RateParam.TimeCol=TimBinFrame+0.5;
                        RateParam.PathSave=[TempResultFolder 'Speed'];
                        RateHist_GroupPlot(TimBinFrame+0.5,GroupSpeed,GroupMetaColor,RateParam)
                        hold on;
                        % text(-10,0.1,['n = ' num2str(CellSampleN(iCell,iPower))])
                        
                        set(gca,'xlim',[-PSTHparam.PreSLMCal PSTHparam.PostSLMCal],'xtick',[-PSTHparam.PreSLMCal:5:PSTHparam.PostSLMCal])
                        % set(gca,'ylim',[-0.1 10])
                        papersizePX=[0 0 12 8];
                        set(gcf, 'PaperUnits', 'centimeters');
                        set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
    
                        % print(gcf, [TempResultFolder 'Speed.svg'], '-dsvg', '-painters');
                        % print(gcf, [TempResultFolder 'Speed.eps'], '-depsc', '-painters');
                        print(gcf, [TempResultFolder 'Speed.tif'], '-dtiffn', '-painters');
                        % 
                        % 
    
                clear GroupParamNet;
                
                GroupParamNet.GroupColor=GroupMetaColor;
                % GroupParamNet.iscell=iscell;
                GroupParamNet.CellN=length(iscell);
                
                % GroupParamNet.SuccTarget=SuccTarget;
                GroupParamNet.ScoreLim=[-0.3 0.3];
                GroupParamNet.ResponseLim=[-0.2 0.2];
                GroupParamNet.ScoreMap=slanCM('wildfire',64);
                GroupParamNet.ScoreLabel='Speed Corr.';
                GroupParamNet.ResponseMap=slanCM('seismic',64);
                GroupParamNet.NodeColor=NodeColor;
                % GroupParamNet.statCellRes=statGroupRes;
                GroupParamNet.crit_pAll=crit_pAll;
                GroupParamNet.SuccAmp=SuccAmp;
                GroupParamNet.xMat=[0.01 0.25 0.25 0.25];
                GroupParamNet.yMat=[0.7 0.7 0.7 0.7];
                GroupParamNet.TargetCellList=TargetCellList;
                GroupParamNet.TargetCellListFunGroup=TargetCellListFunGroup;
                % GroupParamNet.GroupColor=GroupMetaColor;
                
                
                
                close all
                
                
                
                
                crit_pGroup=1;
                GroupParamNet.ScoreLabel='Speed Corr.';
    
                for iFun = 1:length(GroupTargetCellMeta)
                
                        PowerTestAdj=zeros(length(iscell)+1)+NaN;
                        NodeTarget=zeros(length(iscell)+1,1);
                
                        temp1=statGroupRes(iFun).delta;
                        temp2=statGroupRes(iFun).p>crit_pGroup;
                       % temp2=[];
                        temp1(temp2)=NaN;
                        PowerTestAdj(end,1:length(iscell))=temp1;
                
                        NodeTarget(GroupTargetCellMeta{iFun})=temp1([GroupTargetCellMeta{iFun}]);
                        NodeTarget(end)=1;
                        % PowerTestAdj(end,[GroupTargetCell{iFun}])=NaN;
                
                        PowerTestAdj=PowerTestAdj-diag(diag(PowerTestAdj));
                        PowerTestNode=NodeTarget;
                        GroupParamNet.SLMGroup=iFun;
                        GroupParamNet.ResponseLim=[-0.2 0.2]
                        % GroupParamNet.ResponseLim=[-0.3 0.3];
                        GroupParamNet.GroupTargetCell=GroupTargetCellMeta{iFun};
                        GroupParamNet.TargetCellList=GroupTargetCell;
                        % GroupParamNet.TargetCellListFunGroup=TargetCellListFunGroup;
                        GroupParamNet.statCellRes=statGroupRes(iFun);
                
    
                
                        PSTHstruct=PSTHparam;
                        PSTHstruct.Data=GroupResponse{iFun};
                        PSTHstruct.TimeBinFrame=TimBinFrame;
                       [GgraphOut,Res,r,p,tempFig]=ResSLMFunNetwork_ScoreNodeOrder(PSTHstruct,PowerTestAdj,PowerTestNode,GroupParamNet,rSpeed(:,1,1))
                
                       % [GgraphOut,Res,r,p,tempFig]=ResSLMFunNetwork_ScoreNodeOrder(SLMPosInfo,PSTHstruct,PowerTestAdj,PowerTestNode,GroupParamNet,rSpeed(:,1,1))
                       if isempty(GroupParamNet.GroupTargetCell)
                          GgraphOut.p.MarkerSize(end)=0.1;
                          GgraphOut.p.LineStyle='none';
                       end
                        GgraphOut.p.ArrowSize=8;
                        papersizePX=[0 0 50 18];
                        set(gcf, 'PaperUnits', 'centimeters');
                       set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
                
                       axes(tempFig{1});
                       text(TimBinFrame(end),0,['Trial#' num2str(GroupSampleN(iFun))],'verticalalignment','bottom','HorizontalAlignment','right')
                 
    
                        % saveas(gcf,[TempResultFolder GroupLabel{iFun} 'CellGroupSLMVsSpeed'],'tif');
                        % saveas(gcf,[TempResultFolder GroupLabel{iFun} 'CellGroupSLMVsSpeed.eps'],'epsc');
                        % saveas(gcf,[TempResultFolder GroupLabel{iFun} 'CellGroupSLMVsSpeed'],'fig');
                        % print(gcf, [TempResultFolder GroupMetaName{iFun} 'CellGroupSLMVsSpeed.svg'], '-dsvg', '-painters');
                        % print(gcf, [TempResultFolder GroupMetaName{iFun} 'CellGroupSLMVsSpeed.eps'], '-depsc', '-painters');
                        print(gcf, [TempResultFolder GroupMetaName{iFun} 'CellGroupSLMVsSpeed.tif'], '-dtiffn', '-painters');
                
                            close all
    
                   ResParam.Color=[0.1 0.1 0.1];
                   ResParam.Marker='o';
                   ResParam.MarkerSize=12;
                   ResParam.Rtype='spearman';
                   ResParam.xLim=GroupParamNet.ScoreLim;
                   ResParam.yLim=GroupParamNet.ResponseLim;
                   ResParam.xLabel=GroupParamNet.ScoreLabel;
                   ResParam.yLabel='Response';
                   ResParam.ExcludeColor=[1 0 0];
                   ResParam.PlotExclude=0;
                
                        figure;
                        subplotLU(1,4,1,1,P)
                        
                        if iFun<4
                        I0=setdiff(1:length(iscell),GroupTargetCellAll(GroupTargetCellAll(:,2)==iFun,1));
                        else
                        I0=1:length(iscell);
                        end
                
                        [Res,r,p]=LuPairRegressPlot(squeeze(rSpeed(I0,1,1)),statGroupRes(iFun).delta(I0)',ResParam);
                        set(gca,'ylim',[-0.1 0.2])
                        if iFun<4
                        title(['Exclude ' GroupLabel{iFun}])
                        end            
                        subplotLU(1,4,1,2,P)
                        I0=setdiff(1:length(iscell),GroupTargetCellAll(:,1));
                        [Res,r,p]=LuPairRegressPlot(squeeze(rSpeed(I0,1,1)),statGroupRes(iFun).delta(I0)',ResParam);
                        set(gca,'YTickLabel',[])
                        set(gca,'ylim',[-0.1 0.2])
                        title(['Exclude 3 groups'])
                
                        subplotLU(1,4,1,3,P)
                        I1=intersect(find(statGroupRes(iFun).delta>0),I0);
                        [Res,r,p]=LuPairRegressPlot(squeeze(rSpeed(I1,1,1)),statGroupRes(iFun).delta(I1)',ResParam);
                        ylabel('')
                        set(gca,'YTickLabel',[])
                                set(gca,'ylim',[-0.1 0.2])
                        title(['Exclude 3 groups'])
                
                         subplotLU(1,4,1,4,P)
                        I1=intersect(find(statGroupRes(iFun).delta<0),I0);
                        [Res,r,p]=LuPairRegressPlot(squeeze(rSpeed(I1,1,1)),statGroupRes(iFun).delta(I1)',ResParam);
                        ylabel('')
                        set(gca,'YTickLabel',[])
                        set(gca,'ylim',[-0.1 0.2])
                        title(['Exclude 3 groups'])
                
                
                
                
                        papersizePX=[0 0 32 8];
                        set(gcf, 'PaperUnits', 'centimeters');
                        set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
                        LuFontStandard
                        % print(gcf, [TempResultFolder GroupMetaName{iFun} 'SLMVsSpeedRegress.svg'], '-dsvg', '-painters');
                        % print(gcf, [TempResultFolder GroupMetaName{iFun} 'SLMVsSpeedRegress.eps'], '-depsc', '-painters');
                        print(gcf, [TempResultFolder GroupMetaName{iFun} 'SLMVsSpeedRegress.tif'], '-dtiffn', '-painters');
                % print(gcf, [TempResultFolder 'SingleSLMVsSpeed.png'], '-dpng', '-painters');
                
                end
    
                GroupParamNet.ScoreLabel='Stim Corr.';
                GroupParamNet.ScoreLim=[-0.1 0.1];
    
                for iFun = 1:length(GroupTargetCellMeta)
                
                        PowerTestAdj=zeros(length(iscell)+1)+NaN;
                        NodeTarget=zeros(length(iscell)+1,1);
                
                        temp1=statGroupRes(iFun).delta;
                        temp2=statGroupRes(iFun).p>crit_pGroup;
                       % temp2=[];
                        temp1(temp2)=NaN;
                        PowerTestAdj(end,1:length(iscell))=temp1;
                
                        NodeTarget(GroupTargetCellMeta{iFun})=temp1([GroupTargetCellMeta{iFun}]);
                        NodeTarget(end)=1;
                        % PowerTestAdj(end,[GroupTargetCell{iFun}])=NaN;
                
                        PowerTestAdj=PowerTestAdj-diag(diag(PowerTestAdj));
                        PowerTestNode=NodeTarget;
                        GroupParamNet.SLMGroup=iFun;
                        % GroupParamNet.ResponseLim=[-0.3 0.3];
                        GroupParamNet.GroupTargetCell=GroupTargetCellMeta{iFun};
                        GroupParamNet.TargetCellList=GroupTargetCell;
                        % GroupParamNet.TargetCellListFunGroup=TargetCellListFunGroup;
                        GroupParamNet.statCellRes=statGroupRes(iFun);
                
                
                        PSTHstruct=PSTHparam;
                        PSTHstruct.Data=GroupResponse{iFun};
                        PSTHstruct.TimeBinFrame=TimBinFrame;
  
                       [GgraphOut,Res,r,p,tempFig]=ResSLMFunNetwork_ScoreNodeOrder(PSTHstruct,PowerTestAdj,PowerTestNode,GroupParamNet,rStim(:,1,1))

                       % [GgraphOut,Res,r,p,tempFig]=ResSLMFunNetwork_ScoreNodeOrder(SLMPosInfo,PSTHstruct,PowerTestAdj,PowerTestNode,GroupParamNet,rStim(:,1,1))
                       if isempty(GroupParamNet.GroupTargetCell)
                          GgraphOut.p.MarkerSize(end)=0.1;
                          GgraphOut.p.LineStyle='none';
                       end
                        GgraphOut.p.ArrowSize=8;
                        papersizePX=[0 0 50 18];
                        set(gcf, 'PaperUnits', 'centimeters');
                       set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
                        axes(tempFig{1});
                       text(TimBinFrame(end),0,['Trial#' num2str(GroupSampleN(iFun))],'verticalalignment','bottom','HorizontalAlignment','right')

                        % saveas(gcf,[TempResultFolder GroupLabel{iFun} 'CellGroupSLMVsSpeed'],'tif');
                        % saveas(gcf,[TempResultFolder GroupLabel{iFun} 'CellGroupSLMVsSpeed.eps'],'epsc');
                        % saveas(gcf,[TempResultFolder GroupLabel{iFun} 'CellGroupSLMVsSpeed'],'fig');
                        % print(gcf, [TempResultFolder GroupMetaName{iFun} 'CellGroupSLMVsStim.svg'], '-dsvg', '-painters');
                        % print(gcf, [TempResultFolder GroupMetaName{iFun} 'CellGroupSLMVsStim.eps'], '-depsc', '-painters');
                        print(gcf, [TempResultFolder GroupMetaName{iFun} 'CellGroupSLMVsStim.tif'], '-dtiffn', '-painters');
                
                
                   ResParam.Color=[0.1 0.1 0.1];
                   ResParam.Marker='o';
                   ResParam.MarkerSize=12;
                   ResParam.Rtype='spearman';
                   ResParam.xLim=GroupParamNet.ScoreLim;
                   ResParam.yLim=GroupParamNet.ResponseLim;
                   ResParam.xLabel=GroupParamNet.ScoreLabel;
                   ResParam.yLabel='Response';
                   ResParam.ExcludeColor=[1 0 0];
                   ResParam.PlotExclude=0;
                
                        figure;
                        subplotLU(1,4,1,1,P)
                        
                        if iFun<4
                        I0=setdiff(1:length(iscell),GroupTargetCellAll(GroupTargetCellAll(:,2)==iFun,1));
                        else
                        I0=1:length(iscell);
                        end
                
                        [Res,r,p]=LuPairRegressPlot(squeeze(rStim(I0,1,1)),statGroupRes(iFun).delta(I0)',ResParam);
                        set(gca,'ylim',[-0.1 0.2])
                        if iFun<4
                        title(['Exclude ' GroupLabel{iFun}])
                        end  
                
                        subplotLU(1,4,1,2,P)
                        I0=setdiff(1:length(iscell),GroupTargetCellAll(:,1));
                        [Res,r,p]=LuPairRegressPlot(squeeze(rStim(I0,1,1)),statGroupRes(iFun).delta(I0)',ResParam);
                        set(gca,'YTickLabel',[])
                        set(gca,'ylim',[-0.1 0.2])
                        title(['Exclude 3 groups'])
                
                        subplotLU(1,4,1,3,P)
                        I1=intersect(find(statGroupRes(iFun).delta>0),I0);
                        [Res,r,p]=LuPairRegressPlot(squeeze(rStim(I1,1,1)),statGroupRes(iFun).delta(I1)',ResParam);
                        ylabel('')
                        set(gca,'YTickLabel',[])
                                set(gca,'ylim',[-0.1 0.2])
                        title(['Exclude 3 groups'])
                
                         subplotLU(1,4,1,4,P)
                        I1=intersect(find(statGroupRes(iFun).delta<0),I0);
                        [Res,r,p]=LuPairRegressPlot(squeeze(rStim(I1,1,1)),statGroupRes(iFun).delta(I1)',ResParam);
                        ylabel('')
                        set(gca,'YTickLabel',[])
                        set(gca,'ylim',[-0.1 0.2])
                        title(['Exclude 3 groups'])
                
                
                
                
                        papersizePX=[0 0 32 8];
                        set(gcf, 'PaperUnits', 'centimeters');
                        set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
                        LuFontStandard
                        % print(gcf, [TempResultFolder GroupMetaName{iFun} 'SLMVsStimRegress.svg'], '-dsvg', '-painters');
                        % print(gcf, [TempResultFolder GroupMetaName{iFun} 'SLMVsStimRegress.eps'], '-depsc', '-painters');
                        print(gcf, [TempResultFolder GroupMetaName{iFun} 'SLMVsStimRegress.tif'], '-dtiffn', '-painters');
                % print(gcf, [TempResultFolder 'SingleSLMVsSpeed.png'], '-dpng', '-painters');
                close all
                end
    
            end
        % end
    end

end
