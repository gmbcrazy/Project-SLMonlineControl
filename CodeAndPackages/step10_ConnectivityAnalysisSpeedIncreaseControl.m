
clear all
WorkFolder='E:\LuSLMOnlineTest\SL0855-Emx1G6CII-AAV9CAMKII\03062025\';
ResultFolder = Get_ExpDataFolder(WorkFolder, 'Step1', {'Step1Meta.mat','.png'})
load([ResultFolder 'Step1Meta.mat'])
% load('C:\Users\zhangl33\Projects\Project-SLMonlineControl\CodeAndPackages\subfun\Color\colorMapPN3.mat','colorMapPN1')

ResultFolder=[ProcessFolder 'Results\'];
mkdir(ResultFolder)
ResultFolder=[ResultFolder 'Step3\'];
mkdir(ResultFolder)
ResultFolderSpeed=[ResultFolder 'SpeedIncrease\'];
mkdir(ResultFolderSpeed)


iData=1;
PSTHparam.TestStepFrame=3;
[AlignedtempNData,AlignedInfoTable,statCellRes,TargetResponse,TargetCellResP,TargetCellResR,CellSampleN]=Aligned_FromSuite2p(NData{iData},TargetCellList,iscell,Suite2pTable,PVpower,PSTHparam);

GroupLabel={'L','S','N'};
nGroup=length(SLMPosInfo.Group);
GroupColor=[247 150 111;239 109 249;121 247 111]/255;


% TargetCellBand=zeros(length(iscell),1);
% TargetCellBand(TargetCellList)=TargetCellListFunGroup;
% 
NodeColor=repmat([0.9 0.9 0.9],length(iscell),1);

PowerZeroColor=[0.5 0.5 0.5];

% AlignedtempNData=AlignedNData{1};
GroupList=[1 2 3];
PowerZero=[0 1];
PowerZeroLabel={'SLM','FakeSLM'};
VolOut=[0 1];
VolOutLabel={'NoWhisk','Whisk'}
AwakeState=[1 2];
AwakeStateLabel={'Awake','Ane'}

GroupMetaName=[GroupLabel {'FakeSLM'}];
GroupMetaColor=[GroupColor;PowerZeroColor];

GroupTargetCellAll=[];
for iFun=1:length(GroupList)
    GroupTargetCellAll=[GroupTargetCellAll;[GroupTargetCell{iFun}(:) zeros(size(GroupTargetCell{iFun}(:)))+iFun]];
end
GroupTargetCellMeta=[GroupTargetCell {[]}];

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

PostPreDiffSpeedTh=[0 0.5 2 10000];

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

                 SpeedChangeStats{iVol}=ErrorBarPlotLU(1:length(xData),xData,[],GroupMetaColor(GroupfromMetaIndex,:),2,1,[ResultFolderSpeed VolOutLabel{iVol} '.txt'],GroupPair,GroupfromMetaIndex);

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

                        print(gcf, [ResultFolderSpeed 'SpeedChange.svg'], '-dsvg', '-painters');
                        print(gcf, [ResultFolderSpeed 'SpeedChange.eps'], '-depsc', '-painters');
                        print(gcf, [ResultFolderSpeed 'SpeedChange.tif'], '-dtiffn', '-painters');
                
                            close all





PostPreDiffSpeedTh=[0 0.5 2 10000];
for iSpeedTh=1:length(PostPreDiffSpeedTh)
    ResultFolder=[ResultFolderSpeed num2str(PostPreDiffSpeedTh(iSpeedTh)) '\'];
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
                               I1=I1(find((postSLMSpeed-preSLMSpeed)>=-PostPreDiffSpeedTh(iSpeedTh)));
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
                               I1=I1(find((postSLMSpeed-preSLMSpeed)>=-PostPreDiffSpeedTh(iSpeedTh)));
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
                        RateParam.PathSave=[TempResultFolder 'Speed']
    
                        RateHist_GroupPlot(TimBinFrame+0.5,GroupSpeed,GroupMetaColor,RateParam)
                        hold on;
                        % text(-10,0.1,['n = ' num2str(CellSampleN(iCell,iPower))])
                        
                        set(gca,'xlim',[-PSTHparam.PreSLMCal PSTHparam.PostSLMCal],'xtick',[-PSTHparam.PreSLMCal:5:PSTHparam.PostSLMCal])
                        % set(gca,'ylim',[-0.1 10])
                        papersizePX=[0 0 12 8];
                        set(gcf, 'PaperUnits', 'centimeters');
                        set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
    
                        print(gcf, [TempResultFolder 'Speed.svg'], '-dsvg', '-painters');
                        print(gcf, [TempResultFolder 'Speed.eps'], '-depsc', '-painters');
                        print(gcf, [TempResultFolder 'Speed.tif'], '-dtiffn', '-painters');
                        % 
                        % 
    
                clear GroupParamNet;
                
                GroupParamNet.GroupColor=GroupMetaColor;
                GroupParamNet.iscell=iscell;
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
                
                       [GgraphOut,Res,r,p,tempFig]=ResSLMFunNetwork_ScoreNodeOrder(SLMPosInfo,PSTHstruct,PowerTestAdj,PowerTestNode,GroupParamNet,rSpeed(:,1,1))
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
                        print(gcf, [TempResultFolder GroupMetaName{iFun} 'CellGroupSLMVsSpeed.svg'], '-dsvg', '-painters');
                        print(gcf, [TempResultFolder GroupMetaName{iFun} 'CellGroupSLMVsSpeed.eps'], '-depsc', '-painters');
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
                        print(gcf, [TempResultFolder GroupMetaName{iFun} 'SLMVsSpeedRegress.svg'], '-dsvg', '-painters');
                        print(gcf, [TempResultFolder GroupMetaName{iFun} 'SLMVsSpeedRegress.eps'], '-depsc', '-painters');
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
                
                       [GgraphOut,Res,r,p,tempFig]=ResSLMFunNetwork_ScoreNodeOrder(SLMPosInfo,PSTHstruct,PowerTestAdj,PowerTestNode,GroupParamNet,rStim(:,1,1))
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
                        print(gcf, [TempResultFolder GroupMetaName{iFun} 'CellGroupSLMVsStim.svg'], '-dsvg', '-painters');
                        print(gcf, [TempResultFolder GroupMetaName{iFun} 'CellGroupSLMVsStim.eps'], '-depsc', '-painters');
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
                        print(gcf, [TempResultFolder GroupMetaName{iFun} 'SLMVsStimRegress.svg'], '-dsvg', '-painters');
                        print(gcf, [TempResultFolder GroupMetaName{iFun} 'SLMVsStimRegress.eps'], '-depsc', '-painters');
                        print(gcf, [TempResultFolder GroupMetaName{iFun} 'SLMVsStimRegress.tif'], '-dtiffn', '-painters');
                % print(gcf, [TempResultFolder 'SingleSLMVsSpeed.png'], '-dpng', '-painters');
                close all
                end
    
            end
        % end
    end

end


    % for iAwake=1:length(AwakeState)
    %     % for iPower=1:length(PowerZero)
    %         for iVol=1:length(VolOut)
    %             TempResultFolder=[ResultFolder AwakeStateLabel{iAwake} VolOutLabel{iVol} '\'];
    %             mkdir(TempResultFolder);
    %             clear GroupResponse GroupSpeed statGroupRes GroupSampleN
    % 
    %                 for iFun = 1:length(GroupList)
    %                         I1=find(AlignedInfoTable.Group==GroupList(iFun)&AlignedInfoTable.PowerZero==0&AlignedInfoTable.VolOut==VolOut(iVol)&AlignedInfoTable.AwakeState==AwakeState(iAwake));
    %                         % GroupSampleN(iFun)=length(I1);
    % 
    %                            preSLMSpeed=mean(squeeze(AlignedSpeed(1:PSTHparam.PreSLMCal,I1)),1);
    %                            postSLMSpeed=mean(squeeze(AlignedSpeed(PSTHparam.PreSLMCal+[1:TestStepFrame],I1)),1);
    % 
    %                            if iAwake==1
    %                            I1=I1(find((postSLMSpeed-preSLMSpeed)<=0));
    %                            end
    %                            % temp=find((postSLMSpeed-preSLMSpeed)>=0))
    %                            % if length(temp)<=1
    %                            % else
    %                            %    I1=
    %                            % end
    % 
    %                            GroupSampleN(iFun)=length(I1);
    % 
    %                         if length(I1)==1
    %                            GroupResponse{iFun}=squeeze(AlignedtempNData(:,:,I1));
    %                            GroupSpeed{iFun}=AlignedSpeed(:,I1)';
    %                            preSLMdata=squeeze(AlignedtempNData(:,1:PSTHparam.PreSLMCal,I1));
    %                            postSLMdata=squeeze(AlignedtempNData(:,PSTHparam.PreSLMCal+[1:TestStepFrame],I1));
    %                            for jCell=1:length(iscell)
    %                                temp1=preSLMdata(jCell,:,:);
    %                                temp2=postSLMdata(jCell,:,:);
    %                                [~,p(jCell,1),~,t(jCell)]=ttest2(temp2(:),temp1(:));
    %                                temp3(jCell)=mean(temp2(:))-mean(temp1(:));
    %                            end
    %                         elseif length(I1)>1
    %                            GroupResponse{iFun}=nanmean(AlignedtempNData(:,:,I1),3);           
    %                            GroupSpeed{iFun}=AlignedSpeed(:,I1)';
    %                            preSLMdata=squeeze(AlignedtempNData(:,1:PSTHparam.PreSLMCal,I1));
    %                            postSLMdata=squeeze(AlignedtempNData(:,PSTHparam.PreSLMCal+[1:TestStepFrame],I1));
    %                            for jCell=1:length(iscell)
    %                                temp1=preSLMdata(jCell,:,:);
    %                                temp2=postSLMdata(jCell,:,:);
    %                                [~,p(jCell,1),~,t(jCell)]=ttest2(temp2(:),temp1(:));
    %                                temp3(jCell)=mean(temp2(:))-mean(temp1(:));
    %                            end
    %                            statGroupRes(iFun).p=p;
    %                            statGroupRes(iFun).t=t;
    %                            statGroupRes(iFun).delta=temp3;
    % 
    %                            clear p t temp3;
    % 
    % 
    %                         else
    %                            statGroupRes(iFun).p=zeros(size(iscell))+1;
    %                            statGroupRes(iFun).delta=zeros(size(iscell));
    % 
    %                         end
    % 
    %                 end
    % 
    % 
    % 
    % 
    % 
    %                 %%Add Zero power as addtional group
    %                         I1=find(AlignedInfoTable.PowerZero==1&AlignedInfoTable.VolOut==VolOut(iVol)&AlignedInfoTable.AwakeState==AwakeState(iAwake));
    %                         preSLMSpeed=mean(squeeze(AlignedSpeed(1:PSTHparam.PreSLMCal,I1)),1);
    %                         postSLMSpeed=mean(squeeze(AlignedSpeed(PSTHparam.PreSLMCal+[1:TestStepFrame],I1)),1);
    %                         if iAwake==1
    %                            I1=I1(find((postSLMSpeed-preSLMSpeed)<=0));
    %                         end
    % 
    %                         GroupSampleN(iFun+1)=length(I1);
    %                         if length(I1)==1
    %                            GroupResponse{iFun+1}=squeeze(AlignedtempNData(:,:,I1));
    %                            GroupSpeed{iFun+1}=AlignedSpeed(:,I1)';
    %                         elseif length(I1)>1
    %                            GroupResponse{iFun+1}=nanmean(AlignedtempNData(:,:,I1),3);           
    %                            GroupSpeed{iFun+1}=AlignedSpeed(:,I1)';
    %                            preSLMdata=squeeze(AlignedtempNData(:,1:PSTHparam.PreSLMCal,I1));
    %                            postSLMdata=squeeze(AlignedtempNData(:,PSTHparam.PreSLMCal+[1:TestStepFrame],I1));
    %                            for jCell=1:length(iscell)
    %                                temp1=preSLMdata(jCell,:,:);
    %                                temp2=postSLMdata(jCell,:,:);
    %                                [~,p(jCell,1),~,t(jCell)]=ttest2(temp2(:),temp1(:));
    %                                temp3(jCell)=mean(temp2(:))-mean(temp1(:));
    %                            end
    %                            statGroupRes(iFun+1).p=p;
    %                            statGroupRes(iFun+1).t=t;
    %                            statGroupRes(iFun+1).delta=temp3;
    % 
    %                            clear p t temp3;
    % 
    % 
    %                         else
    %                            statGroupRes(iFun).p=zeros(size(iscell))+1;
    %                            statGroupRes(iFun).delta=zeros(size(iscell));
    % 
    %                         end
    %                  %%Add Zero power as addtional group
    % 
    %                     figure;
    %                     % TargetCellList(iCell)
    %                     % subplotLU(length(TargetCellList),length(PVpower),iCell,iPower,P)
    %                     RateParam.TimeCol=TimBinFrame+0.5;
    %                     RateParam.PathSave=[TempResultFolder 'Speed']
    % 
    %                     RateHist_GroupPlot(TimBinFrame+0.5,GroupSpeed,GroupMetaColor,RateParam)
    %                     hold on;
    %                     % text(-10,0.1,['n = ' num2str(CellSampleN(iCell,iPower))])
    % 
    %                     set(gca,'xlim',[-PSTHparam.PreSLMCal PSTHparam.PostSLMCal],'xtick',[-PSTHparam.PreSLMCal:5:PSTHparam.PostSLMCal])
    %                     % set(gca,'ylim',[-0.1 10])
    %                     papersizePX=[0 0 12 8];
    %                     set(gcf, 'PaperUnits', 'centimeters');
    %                     set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
    % 
    %                     print(gcf, [TempResultFolder 'Speed.svg'], '-dsvg', '-painters');
    %                     print(gcf, [TempResultFolder 'Speed.eps'], '-depsc', '-painters');
    %                     print(gcf, [TempResultFolder 'Speed.tif'], '-dtiffn', '-painters');
    %                     % 
    %                     % 
    % 
    %             clear GroupParamNet;
    % 
    %             GroupParamNet.GroupColor=GroupMetaColor;
    %             GroupParamNet.iscell=iscell;
    %             % GroupParamNet.SuccTarget=SuccTarget;
    %             GroupParamNet.ScoreLim=[-0.3 0.3];
    %             GroupParamNet.ResponseLim=[-0.2 0.2];
    %             GroupParamNet.ScoreMap=slanCM('wildfire',64);
    %             GroupParamNet.ScoreLabel='Speed Corr.';
    %             GroupParamNet.ResponseMap=slanCM('seismic',64);
    %             GroupParamNet.NodeColor=NodeColor;
    %             % GroupParamNet.statCellRes=statGroupRes;
    %             GroupParamNet.crit_pAll=crit_pAll;
    %             GroupParamNet.SuccAmp=SuccAmp;
    %             GroupParamNet.xMat=[0.01 0.25 0.25 0.25];
    %             GroupParamNet.yMat=[0.7 0.7 0.7 0.7];
    %             GroupParamNet.TargetCellList=TargetCellList;
    %             GroupParamNet.TargetCellListFunGroup=TargetCellListFunGroup;
    %             % GroupParamNet.GroupColor=GroupMetaColor;
    % 
    % 
    % 
    %             close all
    % 
    % 
    % 
    % 
    %             crit_pGroup=1;
    %             GroupParamNet.ScoreLabel='Speed Corr.';
    % 
    %             for iFun = 1:length(GroupTargetCellMeta)
    % 
    %                     PowerTestAdj=zeros(length(iscell)+1)+NaN;
    %                     NodeTarget=zeros(length(iscell)+1,1);
    % 
    %                     temp1=statGroupRes(iFun).delta;
    %                     temp2=statGroupRes(iFun).p>crit_pGroup;
    %                    % temp2=[];
    %                     temp1(temp2)=NaN;
    %                     PowerTestAdj(end,1:length(iscell))=temp1;
    % 
    %                     NodeTarget(GroupTargetCellMeta{iFun})=temp1([GroupTargetCellMeta{iFun}]);
    %                     NodeTarget(end)=1;
    %                     % PowerTestAdj(end,[GroupTargetCell{iFun}])=NaN;
    % 
    %                     PowerTestAdj=PowerTestAdj-diag(diag(PowerTestAdj));
    %                     PowerTestNode=NodeTarget;
    %                     GroupParamNet.SLMGroup=iFun;
    %                     GroupParamNet.ResponseLim=[-0.2 0.2]
    %                     % GroupParamNet.ResponseLim=[-0.3 0.3];
    %                     GroupParamNet.GroupTargetCell=GroupTargetCellMeta{iFun};
    %                     GroupParamNet.TargetCellList=GroupTargetCell;
    %                     % GroupParamNet.TargetCellListFunGroup=TargetCellListFunGroup;
    %                     GroupParamNet.statCellRes=statGroupRes(iFun);
    % 
    % 
    % 
    %                     PSTHstruct=PSTHparam;
    %                     PSTHstruct.Data=GroupResponse{iFun};
    %                     PSTHstruct.TimeBinFrame=TimBinFrame;
    % 
    %                    [GgraphOut,Res,r,p,tempFig]=ResSLMFunNetwork_ScoreNodeOrder(SLMPosInfo,PSTHstruct,PowerTestAdj,PowerTestNode,GroupParamNet,rSpeed(:,1,1))
    %                    if isempty(GroupParamNet.GroupTargetCell)
    %                       GgraphOut.p.MarkerSize(end)=0.1;
    %                       GgraphOut.p.LineStyle='none';
    %                    end
    %                     GgraphOut.p.ArrowSize=8;
    %                     papersizePX=[0 0 50 18];
    %                     set(gcf, 'PaperUnits', 'centimeters');
    %                    set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
    % 
    %                    axes(tempFig{1});
    %                    text(TimBinFrame(end),0,['Trial#' num2str(GroupSampleN(iFun))],'verticalalignment','bottom','HorizontalAlignment','right')
    % 
    % 
    %                     % saveas(gcf,[TempResultFolder GroupLabel{iFun} 'CellGroupSLMVsSpeed'],'tif');
    %                     % saveas(gcf,[TempResultFolder GroupLabel{iFun} 'CellGroupSLMVsSpeed.eps'],'epsc');
    %                     % saveas(gcf,[TempResultFolder GroupLabel{iFun} 'CellGroupSLMVsSpeed'],'fig');
    %                     print(gcf, [TempResultFolder GroupMetaName{iFun} 'CellGroupSLMVsSpeed.svg'], '-dsvg', '-painters');
    %                     print(gcf, [TempResultFolder GroupMetaName{iFun} 'CellGroupSLMVsSpeed.eps'], '-depsc', '-painters');
    %                     print(gcf, [TempResultFolder GroupMetaName{iFun} 'CellGroupSLMVsSpeed.tif'], '-dtiffn', '-painters');
    % 
    %                         close all
    % 
    %                ResParam.Color=[0.1 0.1 0.1];
    %                ResParam.Marker='o';
    %                ResParam.MarkerSize=12;
    %                ResParam.Rtype='spearman';
    %                ResParam.xLim=GroupParamNet.ScoreLim;
    %                ResParam.yLim=GroupParamNet.ResponseLim;
    %                ResParam.xLabel=GroupParamNet.ScoreLabel;
    %                ResParam.yLabel='Response';
    %                ResParam.ExcludeColor=[1 0 0];
    %                ResParam.PlotExclude=0;
    % 
    %                     figure;
    %                     subplotLU(1,4,1,1,P)
    % 
    %                     if iFun<4
    %                     I0=setdiff(1:length(iscell),GroupTargetCellAll(GroupTargetCellAll(:,2)==iFun,1));
    %                     else
    %                     I0=1:length(iscell);
    %                     end
    % 
    %                     [Res,r,p]=LuPairRegressPlot(squeeze(rSpeed(I0,1,1)),statGroupRes(iFun).delta(I0)',ResParam);
    %                     set(gca,'ylim',[-0.1 0.2])
    %                     if iFun<4
    %                     title(['Exclude ' GroupLabel{iFun}])
    %                     end            
    %                     subplotLU(1,4,1,2,P)
    %                     I0=setdiff(1:length(iscell),GroupTargetCellAll(:,1));
    %                     [Res,r,p]=LuPairRegressPlot(squeeze(rSpeed(I0,1,1)),statGroupRes(iFun).delta(I0)',ResParam);
    %                     set(gca,'YTickLabel',[])
    %                     set(gca,'ylim',[-0.1 0.2])
    %                     title(['Exclude 3 groups'])
    % 
    %                     subplotLU(1,4,1,3,P)
    %                     I1=intersect(find(statGroupRes(iFun).delta>0),I0);
    %                     [Res,r,p]=LuPairRegressPlot(squeeze(rSpeed(I1,1,1)),statGroupRes(iFun).delta(I1)',ResParam);
    %                     ylabel('')
    %                     set(gca,'YTickLabel',[])
    %                             set(gca,'ylim',[-0.1 0.2])
    %                     title(['Exclude 3 groups'])
    % 
    %                      subplotLU(1,4,1,4,P)
    %                     I1=intersect(find(statGroupRes(iFun).delta<0),I0);
    %                     [Res,r,p]=LuPairRegressPlot(squeeze(rSpeed(I1,1,1)),statGroupRes(iFun).delta(I1)',ResParam);
    %                     ylabel('')
    %                     set(gca,'YTickLabel',[])
    %                     set(gca,'ylim',[-0.1 0.2])
    %                     title(['Exclude 3 groups'])
    % 
    % 
    % 
    % 
    %                     papersizePX=[0 0 32 8];
    %                     set(gcf, 'PaperUnits', 'centimeters');
    %                     set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
    %                     LuFontStandard
    %                     print(gcf, [TempResultFolder GroupMetaName{iFun} 'SLMVsSpeedRegress.svg'], '-dsvg', '-painters');
    %                     print(gcf, [TempResultFolder GroupMetaName{iFun} 'SLMVsSpeedRegress.eps'], '-depsc', '-painters');
    %                     print(gcf, [TempResultFolder GroupMetaName{iFun} 'SLMVsSpeedRegress.tif'], '-dtiffn', '-painters');
    %             % print(gcf, [TempResultFolder 'SingleSLMVsSpeed.png'], '-dpng', '-painters');
    % 
    %             end
    % 
    %             GroupParamNet.ScoreLabel='Stim Corr.';
    %             GroupParamNet.ScoreLim=[-0.1 0.1];
    % 
    %             for iFun = 1:length(GroupTargetCellMeta)
    % 
    %                     PowerTestAdj=zeros(length(iscell)+1)+NaN;
    %                     NodeTarget=zeros(length(iscell)+1,1);
    % 
    %                     temp1=statGroupRes(iFun).delta;
    %                     temp2=statGroupRes(iFun).p>crit_pGroup;
    %                    % temp2=[];
    %                     temp1(temp2)=NaN;
    %                     PowerTestAdj(end,1:length(iscell))=temp1;
    % 
    %                     NodeTarget(GroupTargetCellMeta{iFun})=temp1([GroupTargetCellMeta{iFun}]);
    %                     NodeTarget(end)=1;
    %                     % PowerTestAdj(end,[GroupTargetCell{iFun}])=NaN;
    % 
    %                     PowerTestAdj=PowerTestAdj-diag(diag(PowerTestAdj));
    %                     PowerTestNode=NodeTarget;
    %                     GroupParamNet.SLMGroup=iFun;
    %                     % GroupParamNet.ResponseLim=[-0.3 0.3];
    %                     GroupParamNet.GroupTargetCell=GroupTargetCellMeta{iFun};
    %                     GroupParamNet.TargetCellList=GroupTargetCell;
    %                     % GroupParamNet.TargetCellListFunGroup=TargetCellListFunGroup;
    %                     GroupParamNet.statCellRes=statGroupRes(iFun);
    % 
    % 
    %                     PSTHstruct=PSTHparam;
    %                     PSTHstruct.Data=GroupResponse{iFun};
    %                     PSTHstruct.TimeBinFrame=TimBinFrame;
    % 
    %                    [GgraphOut,Res,r,p]=ResSLMFunNetwork_ScoreNodeOrder(SLMPosInfo,PSTHstruct,PowerTestAdj,PowerTestNode,GroupParamNet,rStim(:,1,1))
    %                    if isempty(GroupParamNet.GroupTargetCell)
    %                       GgraphOut.p.MarkerSize(end)=0.1;
    %                       GgraphOut.p.LineStyle='none';
    %                    end
    %                     GgraphOut.p.ArrowSize=8;
    %                     papersizePX=[0 0 50 18];
    %                     set(gcf, 'PaperUnits', 'centimeters');
    %                    set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
    % 
    %                     % saveas(gcf,[TempResultFolder GroupLabel{iFun} 'CellGroupSLMVsSpeed'],'tif');
    %                     % saveas(gcf,[TempResultFolder GroupLabel{iFun} 'CellGroupSLMVsSpeed.eps'],'epsc');
    %                     % saveas(gcf,[TempResultFolder GroupLabel{iFun} 'CellGroupSLMVsSpeed'],'fig');
    %                     print(gcf, [TempResultFolder GroupMetaName{iFun} 'CellGroupSLMVsStim.svg'], '-dsvg', '-painters');
    %                     print(gcf, [TempResultFolder GroupMetaName{iFun} 'CellGroupSLMVsStim.eps'], '-depsc', '-painters');
    %                     print(gcf, [TempResultFolder GroupMetaName{iFun} 'CellGroupSLMVsStim.tif'], '-dtiffn', '-painters');
    % 
    % 
    %                ResParam.Color=[0.1 0.1 0.1];
    %                ResParam.Marker='o';
    %                ResParam.MarkerSize=12;
    %                ResParam.Rtype='spearman';
    %                ResParam.xLim=GroupParamNet.ScoreLim;
    %                ResParam.yLim=GroupParamNet.ResponseLim;
    %                ResParam.xLabel=GroupParamNet.ScoreLabel;
    %                ResParam.yLabel='Response';
    %                ResParam.ExcludeColor=[1 0 0];
    %                ResParam.PlotExclude=0;
    % 
    %                     figure;
    %                     subplotLU(1,4,1,1,P)
    % 
    %                     if iFun<4
    %                     I0=setdiff(1:length(iscell),GroupTargetCellAll(GroupTargetCellAll(:,2)==iFun,1));
    %                     else
    %                     I0=1:length(iscell);
    %                     end
    % 
    %                     [Res,r,p]=LuPairRegressPlot(squeeze(rStim(I0,1,1)),statGroupRes(iFun).delta(I0)',ResParam);
    %                     set(gca,'ylim',[-0.1 0.2])
    %                     if iFun<4
    %                     title(['Exclude ' GroupLabel{iFun}])
    %                     end  
    % 
    %                     subplotLU(1,4,1,2,P)
    %                     I0=setdiff(1:length(iscell),GroupTargetCellAll(:,1));
    %                     [Res,r,p]=LuPairRegressPlot(squeeze(rStim(I0,1,1)),statGroupRes(iFun).delta(I0)',ResParam);
    %                     set(gca,'YTickLabel',[])
    %                     set(gca,'ylim',[-0.1 0.2])
    %                     title(['Exclude 3 groups'])
    % 
    %                     subplotLU(1,4,1,3,P)
    %                     I1=intersect(find(statGroupRes(iFun).delta>0),I0);
    %                     [Res,r,p]=LuPairRegressPlot(squeeze(rStim(I1,1,1)),statGroupRes(iFun).delta(I1)',ResParam);
    %                     ylabel('')
    %                     set(gca,'YTickLabel',[])
    %                             set(gca,'ylim',[-0.1 0.2])
    %                     title(['Exclude 3 groups'])
    % 
    %                      subplotLU(1,4,1,4,P)
    %                     I1=intersect(find(statGroupRes(iFun).delta<0),I0);
    %                     [Res,r,p]=LuPairRegressPlot(squeeze(rStim(I1,1,1)),statGroupRes(iFun).delta(I1)',ResParam);
    %                     ylabel('')
    %                     set(gca,'YTickLabel',[])
    %                     set(gca,'ylim',[-0.1 0.2])
    %                     title(['Exclude 3 groups'])
    % 
    % 
    % 
    % 
    %                     papersizePX=[0 0 32 8];
    %                     set(gcf, 'PaperUnits', 'centimeters');
    %                     set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
    %                     LuFontStandard
    %                     print(gcf, [TempResultFolder GroupMetaName{iFun} 'SLMVsStimRegress.svg'], '-dsvg', '-painters');
    %                     print(gcf, [TempResultFolder GroupMetaName{iFun} 'SLMVsStimRegress.eps'], '-depsc', '-painters');
    %                     print(gcf, [TempResultFolder GroupMetaName{iFun} 'SLMVsStimRegress.tif'], '-dtiffn', '-painters');
    %             % print(gcf, [TempResultFolder 'SingleSLMVsSpeed.png'], '-dpng', '-painters');
    %             close all
    %             end
    % 
    %         end
    %     % end
    % end










































% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % 
% 
% Param.PlotType=3
% Param.statisP=0;
% Param.LegendShow=0;
% Param.Legend=[]
% 
% 
% 
% 
% 
% 
% clear GroupParamNet;
% 
% GroupParamNet.GroupColor=GroupColor;
% GroupParamNet.iscell=iscell;
% % GroupParamNet.SuccTarget=SuccTarget;
% GroupParamNet.ScoreLim=[-0.3 0.3];
% GroupParamNet.ResponseLim=[-0.2 0.2];
% GroupParamNet.ScoreMap=slanCM('wildfire',64);
% GroupParamNet.ScoreLabel='Speed Corr.';
% GroupParamNet.ResponseMap=slanCM('seismic',64);
% GroupParamNet.NodeColor=NodeColor;
% % GroupParamNet.statCellRes=statGroupRes;
% GroupParamNet.crit_pAll=crit_pAll;
% GroupParamNet.SuccAmp=SuccAmp;
% GroupParamNet.xMat=[0.01 0.25 0.25 0.25];
% GroupParamNet.yMat=[0.7 0.7 0.7 0.7];
% GroupParamNet.TargetCellList=TargetCellList;
% GroupParamNet.TargetCellListFunGroup=TargetCellListFunGroup;
% 
% close all
% 
% GroupTargetCellAll=[];
% for iFun=1:length(GroupList)
%     GroupTargetCellAll=[GroupTargetCellAll;[GroupTargetCell{iFun}(:) zeros(size(GroupTargetCell{iFun}(:)))+iFun]];
% end
% 
% 
% crit_pGroup=1;
% GroupParamNet.ScoreLabel='Speed Corr.';
% 
% for iFun = 1:length(GroupTargetCellMeta)
% 
%         PowerTestAdj=zeros(length(iscell)+1)+NaN;
%         NodeTarget=zeros(length(iscell)+1,1);
% 
%         temp1=statGroupRes(iFun).delta;
%         temp2=statGroupRes(iFun).p>crit_pGroup;
%        % temp2=[];
%         temp1(temp2)=NaN;
%         PowerTestAdj(end,1:length(iscell))=temp1;
% 
%         NodeTarget(GroupTargetCellMeta{iFun})=temp1([GroupTargetCellMeta{iFun}]);
%         NodeTarget(end)=1;
%         % PowerTestAdj(end,[GroupTargetCell{iFun}])=NaN;
% 
%         PowerTestAdj=PowerTestAdj-diag(diag(PowerTestAdj));
%         PowerTestNode=NodeTarget;
%         GroupParamNet.SLMGroup=iFun;
%         GroupParamNet.ScoreLim=[-0.4 0.4];
%         GroupParamNet.ResponseLim=[-0.2 0.2]
%         % GroupParamNet.ResponseLim=[-0.3 0.3];
%         GroupParamNet.GroupTargetCell=GroupTargetCellMeta{iFun};
%         GroupParamNet.TargetCellList=GroupTargetCell;
%         % GroupParamNet.TargetCellListFunGroup=TargetCellListFunGroup;
%         GroupParamNet.statCellRes=statGroupRes(iFun);
% 
% 
%         PSTHstruct=PSTHparam;
%         PSTHstruct.Data=GroupResponse{iFun};
%         PSTHstruct.TimeBinFrame=TimBinFrame;
% 
%        [GgraphOut,Res,r,p]=ResSLMFunNetwork_ScoreNodeOrder(SLMPosInfo,PSTHstruct,PowerTestAdj,PowerTestNode,GroupParamNet,rSpeed(:,1,1))
%         GgraphOut.p.ArrowSize=8;
%         papersizePX=[0 0 50 18];
%         set(gcf, 'PaperUnits', 'centimeters');
%        set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
% 
%         % saveas(gcf,[TempResultFolder GroupLabel{iFun} 'CellGroupSLMVsSpeed'],'tif');
%         % saveas(gcf,[TempResultFolder GroupLabel{iFun} 'CellGroupSLMVsSpeed.eps'],'epsc');
%         % saveas(gcf,[TempResultFolder GroupLabel{iFun} 'CellGroupSLMVsSpeed'],'fig');
%         print(gcf, [TempResultFolder GroupLabel{iFun} 'CellGroupSLMVsSpeed.svg'], '-dsvg', '-painters');
%         print(gcf, [TempResultFolder GroupLabel{iFun} 'CellGroupSLMVsSpeed.eps'], '-depsc', '-painters');
%         print(gcf, [TempResultFolder GroupLabel{iFun} 'CellGroupSLMVsSpeed.tif'], '-dtiffn', '-painters');
% 
% 
%    ResParam.Color=[0.1 0.1 0.1];
%    ResParam.Marker='o';
%    ResParam.MarkerSize=12;
%    ResParam.Rtype='spearman';
%    ResParam.xLim=GroupParamNet.ScoreLim;
%    ResParam.yLim=GroupParamNet.ResponseLim;
%    ResParam.xLabel=GroupParamNet.ScoreLabel;
%    ResParam.yLabel='Response';
%    ResParam.ExcludeColor=[1 0 0];
%    ResParam.PlotExclude=0;
% 
%         figure;
%         subplotLU(1,4,1,1,P)
% 
%         if iFun<4
%         I0=setdiff(1:length(iscell),GroupTargetCellAll(GroupTargetCellAll(:,2)==iFun,1));
%         else
%         I0=1:length(iscell);
%         end
% 
%         [Res,r,p]=LuPairRegressPlot(squeeze(rSpeed(I0,1,1)),statGroupRes(iFun).delta(I0)',ResParam);
%         set(gca,'ylim',[-0.1 0.2])
%         if iFun<4
%         title(['Exclude ' GroupLabel{iFun}])
%         end
%         subplotLU(1,4,1,2,P)
%         I0=setdiff(1:length(iscell),GroupTargetCellAll(:,1));
%         [Res,r,p]=LuPairRegressPlot(squeeze(rSpeed(I0,1,1)),statGroupRes(iFun).delta(I0)',ResParam);
%         set(gca,'YTickLabel',[])
%         set(gca,'ylim',[-0.1 0.2])
%         title(['Exclude 3 groups'])
% 
%         subplotLU(1,4,1,3,P)
%         I1=intersect(find(statGroupRes(iFun).delta>0),I0);
%         [Res,r,p]=LuPairRegressPlot(squeeze(rSpeed(I1,1,1)),statGroupRes(iFun).delta(I1)',ResParam);
%         ylabel('')
%         set(gca,'YTickLabel',[])
%                 set(gca,'ylim',[-0.1 0.2])
%         title(['Exclude 3 groups'])
% 
%          subplotLU(1,4,1,4,P)
%         I1=intersect(find(statGroupRes(iFun).delta<0),I0);
%         [Res,r,p]=LuPairRegressPlot(squeeze(rSpeed(I1,1,1)),statGroupRes(iFun).delta(I1)',ResParam);
%         ylabel('')
%         set(gca,'YTickLabel',[])
%         set(gca,'ylim',[-0.1 0.2])
%         title(['Exclude 3 groups'])
% 
% 
% 
% 
%         papersizePX=[0 0 32 8];
%         set(gcf, 'PaperUnits', 'centimeters');
%         set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
%         LuFontStandard
%         print(gcf, [TempResultFolder GroupLabel{iFun} 'SLMVsSpeedRegress.svg'], '-dsvg', '-painters');
%         print(gcf, [TempResultFolder GroupLabel{iFun} 'SLMVsSpeedRegress.eps'], '-depsc', '-painters');
%         print(gcf, [TempResultFolder GroupLabel{iFun} 'SLMVsSpeedRegress.tif'], '-dtiffn', '-painters');
% % print(gcf, [TempResultFolder 'SingleSLMVsSpeed.png'], '-dpng', '-painters');
% 
% end
% 
% 
% 
% close all
% 
% crit_pGroup=1;
% GroupParamNet.ScoreLabel='Stim Corr.';
% 
% for iFun = 1:length(GroupList)
% 
%         PowerTestAdj=zeros(length(iscell)+1)+NaN;
%         NodeTarget=zeros(length(iscell)+1,1);
% 
%         temp1=statGroupRes(iFun).delta;
%         temp2=statGroupRes(iFun).p>crit_pGroup;
%        % temp2=[];
%         temp1(temp2)=NaN;
%         PowerTestAdj(end,1:length(iscell))=temp1;
% 
%         NodeTarget(GroupTargetCell{iFun})=temp1([GroupTargetCell{iFun}]);
%         NodeTarget(end)=1;
%         % PowerTestAdj(end,[GroupTargetCell{iFun}])=NaN;
% 
%         PowerTestAdj=PowerTestAdj-diag(diag(PowerTestAdj));
%         PowerTestNode=NodeTarget;
%         GroupParamNet.SLMGroup=iFun;
%         GroupParamNet.ScoreLim=[-0.1 0.1];
%         GroupParamNet.ResponseLim=[-0.2 0.2]
%         % GroupParamNet.ResponseLim=[-0.3 0.3];
%         GroupParamNet.GroupTargetCell=GroupTargetCell{iFun};
%         GroupParamNet.TargetCellList=GroupTargetCell;
%         % GroupParamNet.TargetCellListFunGroup=TargetCellListFunGroup;
%         GroupParamNet.statCellRes=statGroupRes(iFun);
% 
% 
%         PSTHstruct=PSTHparam;
%         PSTHstruct.Data=GroupResponse{iFun};
%         PSTHstruct.TimeBinFrame=TimBinFrame;
% 
%        [GgraphOut,Res,r,p]=ResSLMFunNetwork_ScoreNodeOrder(SLMPosInfo,PSTHstruct,PowerTestAdj,PowerTestNode,GroupParamNet,rStim(:,1,1))
%         GgraphOut.p.ArrowSize=8;
%         papersizePX=[0 0 50 18];
%         set(gcf, 'PaperUnits', 'centimeters');
%        set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
% 
%         print(gcf, [TempResultFolder GroupLabel{iFun} 'CellGroupSLMVsStim.svg'], '-dsvg', '-painters');
%         print(gcf, [TempResultFolder GroupLabel{iFun} 'CellGroupSLMVsStim.eps'], '-depsc', '-painters');
%         print(gcf, [TempResultFolder GroupLabel{iFun} 'CellGroupSLMVsStim.tif'], '-dtiffn', '-painters');
% 
% 
% 
%    ResParam.Color=[0.1 0.1 0.1];
%    ResParam.Marker='o';
%    ResParam.MarkerSize=12;
%    ResParam.Rtype='spearman';
%    ResParam.xLim=GroupParamNet.ScoreLim;
%    ResParam.yLim=GroupParamNet.ResponseLim;
%    ResParam.xLabel=GroupParamNet.ScoreLabel;
%    ResParam.yLabel='Response';
%    ResParam.ExcludeColor=GroupColor(GroupParamNet.SLMGroup,:);
%    ResParam.PlotExclude=0;
% 
% 
%    P.xLeft=0.1;
% P.yBottom=0.12;
% P.yTop=0.06;
% P.xInt=0.04;
% P.yInt=0.1;
%         figure;
%         subplotLU(1,4,1,1,P)
%         I0=setdiff(1:length(iscell),GroupTargetCellAll(GroupTargetCellAll(:,2)==iFun,1));
%         [Res,r,p]=LuPairRegressPlot(squeeze(rStim(I0,1,1)),statGroupRes(iFun).delta(I0)',ResParam);
%         set(gca,'ylim',[-0.1 0.2])
%         title(['Exclude ' GroupLabel{iFun}])
% 
%         subplotLU(1,4,1,2,P)
%         I0=setdiff(1:length(iscell),GroupTargetCellAll(:,1));
%         [Res,r,p]=LuPairRegressPlot(squeeze(rStim(I0,1,1)),statGroupRes(iFun).delta(I0)',ResParam);
%         set(gca,'YTickLabel',[])
%         set(gca,'ylim',[-0.1 0.2])
%         title(['Exclude 3 groups'])
% 
%         subplotLU(1,4,1,3,P)
%         I1=intersect(find(statGroupRes(iFun).delta>0),I0);
%         [Res,r,p]=LuPairRegressPlot(squeeze(rStim(I1,1,1)),statGroupRes(iFun).delta(I1)',ResParam);
%         ylabel('')
%         set(gca,'YTickLabel',[])
%                 set(gca,'ylim',[-0.1 0.2])
%         title(['Exclude 3 groups'])
% 
%          subplotLU(1,4,1,4,P)
%         I1=intersect(find(statGroupRes(iFun).delta<0),I0);
%         [Res,r,p]=LuPairRegressPlot(squeeze(rStim(I1,1,1)),statGroupRes(iFun).delta(I1)',ResParam);
%         ylabel('')
%         set(gca,'YTickLabel',[])
%         set(gca,'ylim',[-0.1 0.2])
%         title(['Exclude 3 groups'])
% 
%         papersizePX=[0 0 32 8];
%         set(gcf, 'PaperUnits', 'centimeters');
%         set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
%         LuFontStandard
%         print(gcf, [TempResultFolder GroupLabel{iFun} 'SLMVsStimRegress.svg'], '-dsvg', '-painters');
%         print(gcf, [TempResultFolder GroupLabel{iFun} 'SLMVsStimRegress.eps'], '-depsc', '-painters');
%         print(gcf, [TempResultFolder GroupLabel{iFun} 'SLMVsStimRegress.tif'], '-dtiffn', '-painters');
% % print(gcf, [TempResultFolder 'SingleSLMVsSpeed.png'], '-dpng', '-painters');
% 
% end
% 
% 
% close all
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% close all
% crit_pGroup=1;
% GroupParamNet.ScoreLabel='Stim Corr.';
% for iFun = 1:length(GroupList)
% 
%         PowerTestAdj=zeros(length(iscell)+1)+NaN;
%         NodeTarget=zeros(length(iscell)+1,1);
% 
%         temp1=statGroupRes(iFun).delta;
%         temp2=statGroupRes(iFun).p>crit_pGroup;
%        % temp2=[];
%         temp1(temp2)=NaN;
%         PowerTestAdj(end,1:length(iscell))=temp1;
% 
%         NodeTarget(GroupTargetCell{iFun})=temp1([GroupTargetCell{iFun}]);
%         NodeTarget(end)=1;
%         % PowerTestAdj(end,[GroupTargetCell{iFun}])=NaN;
% 
%         PowerTestAdj=PowerTestAdj-diag(diag(PowerTestAdj));
%         PowerTestNode=NodeTarget;
%         GroupParamNet.SLMGroup=iFun;
%         GroupParamNet.ScoreLim=[-0.2 0.2];
%         % GroupParamNet.ResponseLim=[-0.3 0.3];
%         GroupParamNet.GroupTargetCell=GroupTargetCell{iFun};
%         GroupParamNet.TargetCellList=GroupTargetCell;
%         % GroupParamNet.TargetCellListFunGroup=TargetCellListFunGroup;
%         GroupParamNet.statCellRes=statGroupRes(iFun);
% 
% 
%         PSTHstruct=PSTHparam;
%         PSTHstruct.Data=GroupResponse{iFun};
%         PSTHstruct.TimeBinFrame=TimBinFrame;
% 
%        [GgraphOut,Res,r,p]=ResSLMFunNetwork_ScoreNodeOrder(SLMPosInfo,PSTHstruct,PowerTestAdj,PowerTestNode,GroupParamNet,rStim(:,1,1))
%         GgraphOut.p.ArrowSize=8;
%         papersizePX=[0 0 50 18];
%         set(gcf, 'PaperUnits', 'centimeters');
% set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
% 
%         saveas(gcf,[ResultFolder GroupLabel{iFun} 'CellGroupSLMVsStim'],'tif');
%         saveas(gcf,[ResultFolder GroupLabel{iFun} 'CellGroupSLMVsStim.eps'],'epsc');
%         saveas(gcf,[ResultFolder GroupLabel{iFun} 'CellGroupSLMVsStim'],'fig');
% 
% 
% 
% 
% end
% 
% 
% 
% 
% 
% 
% 
% 
% 
% crit_pGroup=1;
% for iFun = 1:length(GroupList)
% 
%         PowerTestAdj=zeros(length(iscell))+NaN;
%         NodeTarget=zeros(length(iscell),1);
% 
%         TargetC=TargetCellList(iCell);
%         temp1=statGroupRes(iFun).delta;
%         temp2=statGroupRes(iFun).p>crit_pGroup;
%        % temp2=[];
%         temp1(temp2)=NaN;
%         PowerTestAdj(GroupTargetCell{iFun}(1),:)=temp1;
%         NodeTarget(GroupTargetCell{iFun})=temp1(GroupTargetCell{iFun});
% 
%         PowerTestAdj=PowerTestAdj-diag(diag(PowerTestAdj));
%         PowerTestNode=NodeTarget;
% 
%        [GgraphOut,Res,r,p]=ResPosNetwork_ScoreNodeOrder(SLMPosInfo,PowerTestAdj,PowerTestNode,ParamNet,rStim(:,1,2))
%         GgraphOut.p.ArrowSize=8;
%         papersizePX=[0 0 30 9]*2;
% 
% end
% 
% 
% TargetCellList=ParamNet.TargetCellList;
% TargetCellListFunGroup=ParamNet.TargetCellListFunGroup;
% 
% 
% 
% crit_pGroup=1;
% for iFun = 1:length(GroupList)
% 
%         PowerTestAdj=zeros(length(iscell)+1)+NaN;
%         NodeTarget=zeros(length(iscell)+1,1);
% 
%         temp1=statGroupRes(iFun).delta;
%         temp2=statGroupRes(iFun).p>crit_pGroup;
%        % temp2=[];
%         temp1(temp2)=NaN;
%         PowerTestAdj(end,1:length(iscell))=temp1;
% 
%         NodeTarget(GroupTargetCell{iFun})=temp1([GroupTargetCell{iFun}]);
%         NodeTarget(end)=1;
%         % PowerTestAdj(end,[GroupTargetCell{iFun}])=NaN;
% 
%         PowerTestAdj=PowerTestAdj-diag(diag(PowerTestAdj));
%         PowerTestNode=NodeTarget;
%         ParamNet.SLMGroup=iFun;
%         ParamNet.ScoreLim=[-0.4 0.4];
%         ParamNet.ResponseLim=[-0.3 0.3];
%         ParamNet.GroupTargetCell=GroupTargetCell{iFun};
%         PSTHstruct=PSTHparam;
%         PSTHstruct.Data=GroupResponse{iFun};
%         PSTHstruct.TimeBinFrame=TimBinFrame;
% 
%        [GgraphOut,Res,r,p]=ResSLMFunNetwork_ScoreNodeOrder(SLMPosInfo,PSTHstruct,PowerTestAdj,PowerTestNode,ParamNet,rSpeed(:,1,2))
%         GgraphOut.p.ArrowSize=8;
%         papersizePX=[0 0 30 9]*2;
% 
% end
% 
% 
% 
%             imagesc(TimBinFrame+0.5,1:length(iscell),CellResponse{iCell,iPower});
%             colormap(colorMapPN1);
%             set(gca,'clim',[-0.2 0.2]);
%             hold on;
%             plot(TimBinFrame(1),TargetCellList(iCell),'g>');
%             set(gca,'xlim',[-PSTHparam.PreSLMCal PSTHparam.PostSLMCal],'xtick',[-PSTHparam.PreSLMCal:5:PSTHparam.PostSLMCal])
% 
% PointList2=find(SLMPosInfo.SLMTable(:,2)>0)';
% SLMPosInfo.SLMRes.*SLMPosInfo.sampleN>=4;
% [~,PowerI]=ismember(SLMPosInfo.SLMTable(:,2),SLMTestInfo.confSet.UncagingLaserPower)
% PVpower=xmlPower2PVpower(SLMTestInfo.confSet.UncagingLaserPower);
% 
% 
% 
% 
% 
% SampleTh=4;
% PointPower=zeros(size(PointList2));
% for iPoint=1:size(PointList2,1)
%     Point=PointList2(iPoint);
%     if PowerI(Point)>0
%        temp1=min(find(SLMPosInfo.SLMRes(Point,:).*SLMPosInfo.sampleN(Point,:)>=SampleTh))
%        % temp2=SLMPosInfo.SLMRes(Point,:)>0
%        temp2=max(find(SLMPosInfo.sampleN(Point,:)>=SampleTh));
%        if ~isempty(temp1)
%        PointPower(Point)=PVpower(temp1);
%        else
%        PointPower(Point)=PVpower(temp2);
%        end
%     end
% end
% PointList3=find(PointPower>0)';
% 
% 
% 
% 
% 
% (PointList2,:)
% 
% nPlane=length(confSet.ETL);
% if (size(CaData.F,2)==sum(SessFileTable.Suite2pTiffNum)/nPlane)
%    disp('Suite2p Data match Tiff table')
% else
%    disp('Warning, suite2p Data does not match Tiff table')
% end
% 
% 
% 
% 
% % % %%Get Initial spontanous recording motion shifts, and move tiff files folder to final data folder for further processing
% % % [FileGenerateInfo,InitialfileList, InitialfileID] = getExpInfoFiles_NonMat(WorkFolder)
% % % [InitialPixShiftFile, Files] = PixShiftLoad(WorkFolder);
% % % copyInitialRecordedFolders(InitialfileList, DataFolder,WorkFolder);
% % % 
% % % load([ConfigFolder '\PreGenerateTseriesMultiZ\SpontBeh5T_Z11Frame550.mat','TSeriesBrukerTBL']);
% % % TSeriesBrukerTBL1=TSeriesBrukerTBL;
% % % load([ConfigFolder 'PreGenerateTseriesMultiZ\Anesthesia5T_Z11Frame550.mat','TSeriesBrukerTBL']);
% % % TSeriesBrukerTBL2=TSeriesBrukerTBL;
% % % clear TSeriesBrukerTBL
% % % TSeriesBrukerTBL=[TSeriesBrukerTBL1 TSeriesBrukerTBL2];
% 
% 
% 
% 
% %%Exclude sessions with non-correct num. of tiff
% [FileGenerateInfo,fileList, fileIDs,tiffNum] = getExpInfoFiles_NonMat(DataFolder);
% 
% PowerTestTiffNum=confSet.Ziteration*confSet.ZRepetition*length(confSet.ETL);
% GroupFunTiffNum=sum(TSeriesBrukerTBL{1}.Reps)*length(confSet.ETL);
% 
% 
% ValidTiffNum=unique(tiffNum(ismember(fileIDs,InitialfileID)|tiffNum==confSet.Ziteration*confSet.ZRepetition*length(confSet.ETL)|tiffNum==sum(TSeriesBrukerTBL{1}.Reps)*length(confSet.ETL)))
% % [FileGenerateInfo,fileList, fileIDs,tiffNum] = getExpInfoFiles_NonMat(DataFolder);
% ismember(tiffNum,ValidTiffNum);
% Invalidfile=~ismember(tiffNum,ValidTiffNum);
% ExFolder=[DataFolder 'ExcludeFolder\'];
% mkdir(ExFolder);
% copyInitialRecordedFolders(fileList(Invalidfile), ExFolder, DataFolder);
% DelFolders(fileList(Invalidfile), DataFolder) 
% 
% %%Exclude sessions with non-correct num. of tiff
% % [FileGenerateInfo,fileList, fileIDs,tiffNum] = getExpInfoFiles_NonMat(DataFolder);
% [MatFile, MatExp] = ExtractExp_FromMat(DataFolder);
% [FileGenerateInfo,fileList, fileIDs,tiffNum] = getExpInfoFiles_NonMat(DataFolder);
% 
% 
% NewData = table(round(InitialfileID(:)), repmat(InitialPixShiftFile, length(InitialfileID), 1), ...
%     'VariableNames', {'FileID', 'motionMed'});
% 
% MatFile = [NewData;MatFile]; 
% 
% %%Check if there is tiff folder where no experimental .mat files is record (due to wrong deletion when recording)
% [~,I1]=setdiff(fileIDs,MatFile.FileID);
% 
% MissingPixShiftFile=[];
% MissingFileID=[];
% if ~isempty(I1)
%     for i=length(I1)
%         [PixShiftFile(i), ~] = PixShiftLoad([DataFolder fileList{I1(i)}]);
%     end
%     MissingFileID=round(fileIDs(I1));
%     NewData = table(MissingFileID(:), PixShiftFile(:), ...
%     'VariableNames', {'FileID', 'motionMed'});
% 
%     MatFile = [NewData;MatFile]; 
% end
% 
% 
% Invalidfile=ismember(fileIDs,MatFile.FileID(MatFile.motionMed>motionTh)');
% copyInitialRecordedFolders(fileList(Invalidfile), ExFolder, DataFolder);
% DelFolders(fileList(Invalidfile), DataFolder);
% 
% 
% [FileGenerateInfo,fileList, fileIDs,tiffNum] = getExpInfoFiles_NonMat(DataFolder);
% PowerTestfileI=find(tiffNum==PowerTestTiffNum)
% GroupFunfileI=find(tiffNum==GroupFunTiffNum)
% 
% PowerTestSLMtbl=[];
% PowerTestSLMNum=[];
% for iFile=1:length(PowerTestfileI)
%     % [OutTBLTemp,XYTrials,PowerWeight]=MPSeqFolder_GroupTargets(Folder,SLMTestInfo)
%     OutTBLTemp=MPSeqFolder_1TargetXNon([DataFolder fileList{PowerTestfileI(iFile)} '\'],[confSet.SLM_Pixels_Y;confSet.SLM_Pixels_X],SLMTestInfo.Pos3Dneed);
%     OutTBLTemp.FileID = fileIDs(PowerTestfileI(iFile)) * ones(size(OutTBLTemp, 1), 1);
%     PowerTestSLMtbl=[PowerTestSLMtbl;OutTBLTemp];
%     PowerTestSLMNum(iFile,:)=[fileIDs(PowerTestfileI(iFile)) size(OutTBLTemp,1)];
% end
% 
% 
% Pos3DFun=SLMPosInfo.FinalPos3D;
% % FunScore=SLMPosInfo.FinalFunScore;
% Group=SLMPosInfo.Group;
% for iFun=1:length(Group)
%     Pos3DGroup{iFun}=Pos3DFun(Group(iFun).Indices,:);
% end
% 
% 
% FunSLMtbl=[];
% FunSLMNum=[];
% for iFile=1:length(GroupFunfileI)
%     % [OutTBLTemp,XYTrials,PowerWeight]=MPSeqFolder_GroupTargets(Folder,SLMTestInfo)
%     OutTBLTemp=MPSeqFolder_GroupTargets([DataFolder fileList{GroupFunfileI(iFile)} '\'],[confSet.SLM_Pixels_Y;confSet.SLM_Pixels_X],Pos3DGroup);
%     OutTBLTemp.FileID = fileIDs(GroupFunfileI(iFile)) * ones(size(OutTBLTemp, 1), 1);
%     FunSLMtbl=[FunSLMtbl;OutTBLTemp];
%     FunSLMNum(iFile)=size(OutTBLTemp,1);
% end
% 
% 
% InvalidI1=PowerTestfileI(PowerTestSLMNum(:,2)<(confSet.Ziteration-1));
% InvalidI2=PowerTestfileI(FunSLMNum(:,2)<sum(TSeriesBrukerTBL{1}.SynMP));
% InvalidI=union(InvalidI1,InvalidI2);
% copyInitialRecordedFolders(fileList(InvalidI), ExFolder, DataFolder);
% DelFolders(fileList(InvalidI), DataFolder);
% 
% 
% FunSLMtbl = MatchOutTBLAll_TSeriesBruker(FunSLMtbl, TSeriesBrukerTBL);
% 
% T = outerjoin(PowerTestSLMtbl, FunSLMtbl, 'MergeKeys', true);
% [~,I1]=sort(T.FileID);
% T=T(I1,:);
% 
% FileInfo=table(fileIDs(:),tiffNum(:), fileList(:),'VariableNames',{'FileID', 'tiffNum','FileKey'});
% 
% FileInfo = innerjoin(FileInfo, MatFile, "Keys", "FileID");
% 
% TT = outerjoin(T, FileInfo, "Keys", "FileID", "Type", "right", "MergeKeys", true);
% 
% 
% 
% OutTBLAll=TT; 
% OutTBLAll.AwakeState(TT.TSeriesInd<=5)=1;    %%The 1st half 5 Tseries is designed for awake state.
% OutTBLAll.AwakeState(TT.TSeriesInd>=6)=2;    %%The 2nd half 5 Tseries is designed for anesia state.
% 
% 
% FileID=unique(OutTBLAll.FileID)
% 
% [~,~, fileIDCurrent,tiffNum] = getExpInfoFiles_NonMat(DataFolder);
% if isempty(setdiff(fileIDCurrent,FileID))&&isempty(setdiff(FileID,fileIDCurrent))
%    disp('FileID in tiff folder and OutTBLAll Table match, continue to delete post SLM tiff files for all folders');
%    RemoveFrame=2; %%2 repetitions right together with MarkPoint would be removed.
%     [TiffTable, RemoveList] = RemoveMPsynTiffFolder(DataFolder,RemoveFrame);
% 
%     [FileGenerateInfo,fileList, fileIDs,tiffNum] = getExpInfoFiles_NonMat(DataFolder);
%     PostTiffTable=table(fileIDs(:),tiffNum(:),'VariableNames',{'FileID','Suite2pTiffNum'})
%     Suite2pTable=outerjoin(OutTBLAll, PostTiffTable, "Keys", "FileID", "Type", "right", "MergeKeys", true); 
%     save([DataFolder 'TableForSuite2p.mat'],'Suite2pTable','SLMPosInfo','SLMTestInfo');
% 
%     [~,i1]=unique(Suite2pTable.FileID);
%     SessFileTable=Suite2pTable(i1,:);
%     disp(['Total of ' num2str(sum(SessFileTable.Suite2pTiffNum)) ' tif files required to processed in suite2p']);
% else
%     disp('FileID in tiff folder and OutTBLAll Table do NOT match, check!');
% end
% 
% sponRecording = findAllFoldersKeyWords(WorkFolder, 'TSeries',0);
% 
