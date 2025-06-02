clear all
BatchSavePath='D:\Project1-LocalProcessing\Step1\';
load([BatchSavePath 'FOV.mat'])
Suite2pDataKeywords='awakeRefSpon';

DataSavePath='D:\Project1-LocalProcessing\Step3\';
mkdir(DataSavePath);
DataSavePath=[DataSavePath Suite2pDataKeywords '\'];
mkdir(DataSavePath);
SaveFunCon=[DataSavePath 'FunCon\'];
mkdir(SaveFunCon)



PSTHparam.PreSLMCal = 10; 
PSTHparam.PostSLMCal = 3;
PSTHparam.pTh = 0.05; 
PSTHparam.TestMethod = 'ranksum';
PSTHparam.MPFrameJump = 2;
PSTHparam.TestStepFrame = 3;    %%post-slm frames for Test whether SLM works
PSTHparam.iData = 1;    %%post-slm frames for Test whether SLM works


TimBinFrame = -PSTHparam.PreSLMCal:PSTHparam.PostSLMCal-1;
%% Initial align behaviors with imaging, identify SLM target cells
% for iFOV=1:length(FOVUpdate)
% % iFOV=3;
%     FOVtemp=FOVUpdate(iFOV);
%     suite2pFOVPathLocalTemp=suite2pFOVPathLocal{iFOV};
%     OfflineCombineProcess_OneFOV(FOVtemp, Suite2pDataKeywords, suite2pFOVPathLocalTemp,PSTHparam);
% end
% NodeColor=repmat([0.9 0.9 0.9],NCell,1);

% %% Power test data
% for iFOV=1:length(FOVUpdate)
% % iFOV=3;
%     FOVtemp=FOVUpdate(iFOV);
%     suite2pFOVPathLocalTemp=suite2pFOVPathLocal{iFOV};
%     OfflinePowerTest_OneFOV(FOVtemp, Suite2pDataKeywords, suite2pFOVPathLocalTemp,PSTHparam,SaveFunCon)
% end




%% Group SLM data
SaveFunCon=[DataSavePath 'GroupSLM6Sessions\'];
mkdir(SaveFunCon)
% for iFOV=1:length(FOVUpdate)
%     FOVtemp=FOVUpdate(iFOV);
%     suite2pFOVPathLocalTemp=suite2pFOVPathLocal{iFOV};
%     OfflineSLMAbsSpeedControl_OneFOV(FOVtemp, Suite2pDataKeywords, suite2pFOVPathLocalTemp,PSTHparam,SaveFunCon)
% end
IndexFOVNeed=[1 2 3 4 7 8];
Output=OfflineSLM_ExtractFOVs(FOVUpdate(IndexFOVNeed), Suite2pDataKeywords,suite2pFOVPathLocal(IndexFOVNeed),PSTHparam);


GroupLabel={'L','S','N'};
GroupList=[1 2 3];

for iFun=1:length(GroupList)
    GroupTargetCellTemp=[];
    Cnum=0;
    for iFOV = 1:length(Output.GroupTargetCellMeta)
         GroupTargetCellTemp=[GroupTargetCellTemp;Output.GroupTargetCellMeta{iFOV}{iFun}(:)+Cnum];
         Cnum = Cnum+sum(abs(Output.NeuroPos3DMeta(:,4)-iFOV)<0.1);
    end
    GroupTargetCell{iFun}=GroupTargetCellTemp;
end
GroupTargetCellAll=[];
for iFun=1:length(GroupList)
    GroupTargetCellAll=[GroupTargetCellAll;[GroupTargetCell{iFun}(:) zeros(size(GroupTargetCell{iFun}(:)))+iFun]];
end
GroupTargetCellMeta=[GroupTargetCell {[]}];

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

SaveP1=SaveFunCon


AlignedInfoTable=Output.AlignedInfoTable;
AlignedSpeed=Output.AlignedSpeedMeta;
% AlignedSpeed=Output.AlignedSpeed;
AlignedtempNData=Output.AlignedNData;

PSTHparam.TestStepFrame=3;
TestStepFrame=PSTHparam.TestStepFrame;

PostPreDiffSpeedTh=[1 2 10000];

NCell=size(Output.NeuroPos3DMeta,1);
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


rSpeed=Output.rSpeed;
rStim=Output.rStim;

for iSpeedTh=1:length(PostPreDiffSpeedTh)
    ResultFolder=[SaveP1 num2str(PostPreDiffSpeedTh(iSpeedTh)) '\'];
    mkdir(ResultFolder);

    for iAwake=1:length(AwakeState)
        % for iPower=1:length(PowerZero)
            for iVol=1:length(VolOut)
                TempResultFolder=[ResultFolder AwakeStateLabel{iAwake} VolOutLabel{iVol} '\'];
                mkdir(TempResultFolder);
                clear statGroupRes GroupSampleN GroupSpeed GroupResponse
                % clear GroupResponse GroupSpeed statGroupRes GroupSampleN
  
                statGroupRes(iFun).p=zeros(NCell,1)+1;
                statGroupRes(iFun).delta=zeros(NCell,1);
                GroupResponse{iFun}=zeros(NCell,length(TimBinFrame));

                    for iFun = 1:length(GroupList)
                            I1All=find(AlignedInfoTable.Group==GroupList(iFun)&AlignedInfoTable.PowerZero==0&AlignedInfoTable.VolOut==VolOut(iVol)&AlignedInfoTable.AwakeState==AwakeState(iAwake));
                            % GroupSampleN(iFun)=length(I1);
    
                               preSLMSpeed=mean(squeeze(AlignedSpeed(1:PSTHparam.PreSLMCal,I1All)),1);
                               postSLMSpeed=mean(squeeze(AlignedSpeed(PSTHparam.PreSLMCal+[1:TestStepFrame],I1All)),1);
    
                               if iAwake==1
                               I1All=I1All(find(abs(postSLMSpeed-preSLMSpeed)<=PostPreDiffSpeedTh(iSpeedTh)));
                               end
                               % temp=find((postSLMSpeed-preSLMSpeed)>=0))
                               % if length(temp)<=1
                               % else
                               %    I1=
                               % end
                               GroupSpeed{iFun}=AlignedSpeed(:,I1All)';
                               GroupSampleN(iFun)=length(I1All);
    

                            for iFOV = 1:length(Output.AlignedNData)
                                AlignedInfoTableFOV=Output.AlignedInfoTableFOV{iFOV};
                                I1=find(AlignedInfoTableFOV.Group==GroupList(iFun)&AlignedInfoTableFOV.PowerZero==0&AlignedInfoTableFOV.VolOut==VolOut(iVol)&AlignedInfoTableFOV.AwakeState==AwakeState(iAwake));
                                AlignedtempNData = Output.AlignedNData{iFOV};
                                iFOVneuro=find(Output.NeuroPos3DMeta(:,4)==iFOV);

                                iFOVEvent = find(AlignedInfoTable.iFOV==iFOV);

                                AlignedSpeedFOV=AlignedSpeed(:,iFOVEvent);

                                preSLMSpeed=mean(squeeze(AlignedSpeedFOV(1:PSTHparam.PreSLMCal,I1)),1);
                                postSLMSpeed=mean(squeeze(AlignedSpeedFOV(PSTHparam.PreSLMCal+[1:TestStepFrame],I1)),1);
                                if iAwake==1
                                   I1=I1(find(abs(postSLMSpeed-preSLMSpeed)<=PostPreDiffSpeedTh(iSpeedTh)));
                                end

                                if length(I1)==1
                                   % GroupResponse{iFun}=squeeze(AlignedtempNData(:,:,I1));
                                   % GroupSpeed{iFun}=AlignedSpeed(:,I1)';
                                   GroupResponse{iFun}(iFOVneuro,:)=squeeze(AlignedtempNData(:,:,I1));
                                   preSLMdata=squeeze(AlignedtempNData(:,1:PSTHparam.PreSLMCal,I1));
                                   postSLMdata=squeeze(AlignedtempNData(:,PSTHparam.PreSLMCal+[1:TestStepFrame],I1));
                                   for jCell=1:sum(Output.NeuroPos3DMeta(:,4)==iFOV)
                                       temp1=preSLMdata(jCell,:,:);
                                       temp2=postSLMdata(jCell,:,:);
                                       [~,p(jCell,1),~,t(jCell)]=ttest2(temp2(:),temp1(:));
                                       temp3(jCell)=mean(temp2(:))-mean(temp1(:));
                                   end
                                elseif length(I1)>1
                                   % GroupResponse{iFun}=nanmean(AlignedtempNData(:,:,I1),3);           
                                   % GroupSpeed{iFun}=AlignedSpeed(:,I1)';
                                   GroupResponse{iFun}(iFOVneuro,:)=nanmean(AlignedtempNData(:,:,I1),3);  
                                   preSLMdata=squeeze(AlignedtempNData(:,1:PSTHparam.PreSLMCal,I1));
                                   postSLMdata=squeeze(AlignedtempNData(:,PSTHparam.PreSLMCal+[1:TestStepFrame],I1));
                                   for jCell=1:sum(Output.NeuroPos3DMeta(:,4)==iFOV)
                                       temp1=preSLMdata(jCell,:,:);
                                       temp2=postSLMdata(jCell,:,:);
                                       [~,p(jCell,1),~,t(jCell)]=ttest2(temp2(:),temp1(:));
                                       temp3(jCell)=mean(temp2(:))-mean(temp1(:));
                                   end
                                   statGroupRes(iFun).p(iFOVneuro)=p;
                                   % statGroupRes(iFun).t=t;
                                   statGroupRes(iFun).delta(iFOVneuro)=temp3;
                        
                                   clear p t temp3;
                        
                        
                                else
                                   % statGroupRes(iFun).p=zeros(size(iscell))+1;
                                   % statGroupRes(iFun).delta=zeros(size(iscell));
        
                                end

                            end







                    
                    end
    
    
    
                statGroupRes(iFun+1).p=zeros(NCell,1)+1;
                statGroupRes(iFun+1).delta=zeros(NCell,1);
                GroupResponse{iFun+1}=zeros(NCell,length(TimBinFrame));

    
                    %%Add Zero power as addtional group
                            I1All=find(AlignedInfoTable.PowerZero==1&AlignedInfoTable.VolOut==VolOut(iVol)&AlignedInfoTable.AwakeState==AwakeState(iAwake));
                            preSLMSpeed=mean(squeeze(AlignedSpeed(1:PSTHparam.PreSLMCal,I1All)),1);
                            postSLMSpeed=mean(squeeze(AlignedSpeed(PSTHparam.PreSLMCal+[1:TestStepFrame],I1All)),1);
                            if iAwake==1
                               I1All=I1All(find(abs(postSLMSpeed-preSLMSpeed)<=PostPreDiffSpeedTh(iSpeedTh)));
                            end
                            GroupSpeed{iFun+1}=AlignedSpeed(:,I1All)';

                            GroupSampleN(iFun+1)=length(I1All);

                            for iFOV = 1:length(Output.AlignedNData)
                                    AlignedInfoTableFOV=Output.AlignedInfoTableFOV{iFOV};
                                    I1=find(AlignedInfoTableFOV.PowerZero==1&AlignedInfoTableFOV.VolOut==VolOut(iVol)&AlignedInfoTableFOV.AwakeState==AwakeState(iAwake));
                                    AlignedtempNData = Output.AlignedNData{iFOV};
                                    iFOVneuro=find(Output.NeuroPos3DMeta(:,4)==iFOV);
                                if length(I1)==1
                                   % GroupResponse{iFun+1}=squeeze(AlignedtempNData(:,:,I1));
                                   % GroupSpeed{iFun+1}=AlignedSpeed(:,I1)';
                                   GroupResponse{iFun+1}(iFOVneuro,:)=squeeze(AlignedtempNData(:,:,I1));

                                elseif length(I1)>1
                                   % GroupResponse{iFun+1}=nanmean(AlignedtempNData(:,:,I1),3);    
                                   GroupResponse{iFun+1}(iFOVneuro,:)=nanmean(AlignedtempNData(:,:,I1),3);  
                                   % GroupSpeed{iFun+1}=AlignedSpeed(:,I1)';
                                   preSLMdata=squeeze(AlignedtempNData(:,1:PSTHparam.PreSLMCal,I1));
                                   postSLMdata=squeeze(AlignedtempNData(:,PSTHparam.PreSLMCal+[1:TestStepFrame],I1));
                                   for jCell=1:sum(Output.NeuroPos3DMeta(:,4)==iFOV)
                                       temp1=preSLMdata(jCell,:,:);
                                       temp2=postSLMdata(jCell,:,:);
                                       [~,p(jCell,1),~,t(jCell)]=ttest2(temp2(:),temp1(:));
                                       temp3(jCell)=mean(temp2(:))-mean(temp1(:));
                                   end
                                   statGroupRes(iFun+1).p(iFOVneuro)=p;
                                   % statGroupRes(iFun+1).t=t;
                                   statGroupRes(iFun+1).delta(iFOVneuro)=temp3;
                        
                                   clear p t temp3;
                        
                        
                                else
                                   % statGroupRes(iFun+1).p=zeros(size(iscell))+1;
                                   % statGroupRes(iFun+1).delta=zeros(size(iscell));
        
                                end

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
                % GroupParamNet.SuccTarget=SuccTarget;
                GroupParamNet.ScoreLim=[-0.3 0.3];
                GroupParamNet.ResponseLim=[-0.2 0.2];
                GroupParamNet.ScoreMap=slanCM('wildfire',64);
                GroupParamNet.ScoreLabel='Speed Corr.';
                GroupParamNet.ResponseMap=slanCM('seismic',64);
                GroupParamNet.NodeColor=NodeColor;
                % GroupParamNet.statCellRes=statGroupRes;
                GroupParamNet.crit_pAll=1;
                % GroupParamNet.SuccAmp=SuccAmp;
                GroupParamNet.xMat=[0.01 0.25 0.25 0.25];
                GroupParamNet.yMat=[0.7 0.7 0.7 0.7];
                % GroupParamNet.TargetCellList=TargetCellList;
                % GroupParamNet.TargetCellListFunGroup=TargetCellListFunGroup;
                GroupParamNet.CellN=NCell;
                GroupParamNet.FOVMap=slanCM("Paired",length(FOVUpdate));
                GroupParamNet.FOVLabel='Session';
                GroupParamNet.Data=GroupResponse;
                % GroupParamNet.GroupColor=GroupMetaColor;
                
                
                
                close all
                
                
                
                
                crit_pGroup=1;
                GroupParamNet.ScoreLabel='Speed Corr.';
    
                for iFun = 1:length(GroupTargetCellMeta)
                
                        PowerTestAdj=zeros(NCell+1)+NaN;
                        NodeTarget=zeros(NCell+1,1);
                
                        temp1=statGroupRes(iFun).delta;
                        temp2=statGroupRes(iFun).p>crit_pGroup;
                       % temp2=[];
                        temp1(temp2)=NaN;
                        PowerTestAdj(end,1:NCell)=temp1;
                
                        NodeTarget(GroupTargetCellMeta{iFun})=temp1([GroupTargetCellMeta{iFun}]);
                        NodeTarget(end)=1;
                        % PowerTestAdj(end,[GroupTargetCell{iFun}])=NaN;
                
                        PowerTestAdj=PowerTestAdj-diag(diag(PowerTestAdj));
                        PowerTestNode=NodeTarget;
                        GroupParamNet.SLMGroup=iFun;
                        GroupParamNet.ResponseLim=[-0.2 0.2]
                        % GroupParamNet.ResponseLim=[-0.3 0.3];
                        GroupParamNet.GroupTargetCell=GroupTargetCellMeta{iFun};
                        GroupParamNet.TargetCellList=GroupTargetCellMeta{iFun};

                        % GroupParamNet.TargetCellListFunGroup=TargetCellListFunGroup;
                        GroupParamNet.statCellRes=statGroupRes(iFun);
                        GroupParamNet.xMat=[0.01 0.25 0.25 0.25];
                        GroupParamNet.yMat=[0.7 0.7 0.7 0.7];
                        GroupParamNet.CellN=NCell;
                         % GroupParamNet.Data=GroupResponse{iFun};

                        % GroupParamNet.TargetCellList=GroupTargetCellMeta{iFun};
                        % GroupParamNet.TargetCellListFunGroup=TargetCellListFunGroup;
                        % GroupParamNet.FOVMap=slanCM("Pair");

                
                        PSTHstruct=PSTHparam;

                       PSTHstruct.Data=GroupResponse{iFun};
                       PSTHstruct.TimeBinFrame=TimBinFrame;
                
                       [GgraphOut,Res,r,p,tempFig]=ResSLMFunNetwork_ScoreNodeOrder_MetaFOV(PSTHstruct,PowerTestAdj,PowerTestNode,GroupParamNet,rSpeed(:,1,1),Output.NeuroPos3DMeta(:,4));
                       if isempty(GroupParamNet.GroupTargetCell)
                          GgraphOut.p.MarkerSize(end)=0.1;
                          GgraphOut.p.LineStyle='none';
                       end
                       % GgraphOut.p.NodeColor=NodeColorSorted;
                       % theta=ClockWiseGraph1(GgraphOut.G,GgraphOut.p,radius);


                        GgraphOut.p.ArrowSize=6;
                        papersizePX=[0 0 50 18];
                        set(gcf, 'PaperUnits', 'centimeters');
                       set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
                
                       axes(tempFig{1});
                       text(TimBinFrame(end),0,['Trial#' num2str(GroupSampleN(iFun))],'verticalalignment','bottom','HorizontalAlignment','right')
                 
    
                        % saveas(gcf,[TempResultFolder GroupLabel{iFun} 'CellGroupSLMVsSpeed'],'tif');
                        % saveas(gcf,[TempResultFolder GroupLabel{iFun} 'CellGroupSLMVsSpeed.eps'],'epsc');
                        % saveas(gcf,[TempResultFolder GroupLabel{iFun} 'CellGroupSLMVsSpeed'],'fig');
                        print(gcf, [TempResultFolder GroupMetaName{iFun} 'CellGroupSLMVsSpeed.svg'], '-dsvg', '-painters');
                        % print(gcf, [TempResultFolder GroupMetaName{iFun} 'CellGroupSLMVsSpeed.eps'], '-depsc', '-painters');
                        saveas(gcf, [TempResultFolder GroupMetaName{iFun} 'CellGroupSLMVsSpeed'], 'tif');
                
                        close all
    
                       ResParam.Color=GroupParamNet.FOVMap;
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
                        I0=setdiff(1:NCell,GroupTargetCellAll(GroupTargetCellAll(:,2)==iFun,1));
                        else
                        I0=1:NCell;
                        end
                
                        [Res,r,p]=LuPairRegressPlot_Group(squeeze(rSpeed(I0,1,1)),statGroupRes(iFun).delta(I0)',Output.NeuroPos3DMeta(I0,4),ResParam,statGroupRes(4).delta(I0));
                        set(gca,'ylim',[-0.1 0.2])
                        if iFun<4
                        title(['Exclude ' GroupLabel{iFun}])
                        end            
                        subplotLU(1,4,1,2,P)
                        I0=setdiff(1:NCell,GroupTargetCellAll(:,1));
                        [Res,r,p]=LuPairRegressPlot_Group(squeeze(rSpeed(I0,1,1)),statGroupRes(iFun).delta(I0)',Output.NeuroPos3DMeta(I0,4),ResParam,statGroupRes(4).delta(I0));
                        set(gca,'YTickLabel',[])
                        set(gca,'ylim',[-0.1 0.2])
                        title(['Exclude 3 groups'])
                
                        subplotLU(1,4,1,3,P)
                        I1=intersect(find(statGroupRes(iFun).delta>0),I0);
                        [Res,r,p]=LuPairRegressPlot_Group(squeeze(rSpeed(I1,1,1)),statGroupRes(iFun).delta(I1)',Output.NeuroPos3DMeta(I1,4),ResParam,statGroupRes(4).delta(I1));
                        ylabel('')
                        set(gca,'YTickLabel',[])
                                set(gca,'ylim',[-0.1 0.2])
                        title(['Exclude 3 groups'])
                
                         subplotLU(1,4,1,4,P)
                        I1=intersect(find(statGroupRes(iFun).delta<0),I0);
                        [Res,r,p]=LuPairRegressPlot_Group(squeeze(rSpeed(I1,1,1)),statGroupRes(iFun).delta(I1)',Output.NeuroPos3DMeta(I1,4),ResParam,statGroupRes(4).delta(I1));
                        ylabel('')
                        set(gca,'YTickLabel',[])
                        set(gca,'ylim',[-0.1 0.2])
                        title(['Exclude 3 groups'])
                
                
                        clear p
                
                        papersizePX=[0 0 32 8];
                        set(gcf, 'PaperUnits', 'centimeters');
                        set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
                        LuFontStandard
                        print(gcf, [TempResultFolder GroupMetaName{iFun} 'SLMVsSpeedRegress.svg'], '-dsvg', '-painters');
                        % print(gcf, [TempResultFolder GroupMetaName{iFun} 'SLMVsSpeedRegress.eps'], '-depsc', '-painters');
                        saveas(gcf, [TempResultFolder GroupMetaName{iFun} 'SLMVsSpeedRegress'], 'tif');
                % print(gcf, [TempResultFolder 'SingleSLMVsSpeed.png'], '-dpng', '-painters');
                
                end
    
                GroupParamNet.ScoreLabel='Stim Corr.';
                GroupParamNet.ScoreLim=[-0.1 0.1];
    
                for iFun = 1:length(GroupTargetCellMeta)
                
                        PowerTestAdj=zeros(NCell+1)+NaN;
                        NodeTarget=zeros(NCell+1,1);
                
                        temp1=statGroupRes(iFun).delta;
                        temp2=statGroupRes(iFun).p>crit_pGroup;
                       % temp2=[];
                        temp1(temp2)=NaN;
                        PowerTestAdj(end,1:NCell)=temp1;
                
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
                
                       [GgraphOut,Res,r,p,tempFig]=ResSLMFunNetwork_ScoreNodeOrder_MetaFOV(PSTHstruct,PowerTestAdj,PowerTestNode,GroupParamNet,rStim(:,1,1),Output.NeuroPos3DMeta(:,4));

                       % [GgraphOut,Res,r,p,tempFig]=ResSLMFunNetwork_ScoreNodeOrder(PSTHstruct,PowerTestAdj,PowerTestNode,GroupParamNet,rStim(:,1,1))
                       % if isempty(GroupParamNet.GroupTargetCell)
                       %    GgraphOut.p.MarkerSize(end)=0.1;
                       %    GgraphOut.p.LineStyle='none';
                       % end
                       GgraphOut.p.ArrowSize=6;
                        papersizePX=[0 0 50 18];
                        set(gcf, 'PaperUnits', 'centimeters');
                       set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
                        axes(tempFig{1});
                       text(TimBinFrame(end),0,['Trial#' num2str(GroupSampleN(iFun))],'verticalalignment','bottom','HorizontalAlignment','right')

                        % saveas(gcf,[TempResultFolder GroupLabel{iFun} 'CellGroupSLMVsSpeed'],'tif');
                        % saveas(gcf,[TempResultFolder GroupLabel{iFun} 'CellGroupSLMVsSpeed.eps'],'epsc');
                        % saveas(gcf,[TempResultFolder GroupLabel{iFun} 'CellGroupSLMVsSpeed'],'fig');
                        print(gcf, [TempResultFolder GroupMetaName{iFun} 'CellGroupSLMVsStim.svg'], '-dsvg', '-painters');
                        % print(gcf, [TempResultFolder GroupMetaName{iFun} 'CellGroupSLMVsStim.eps'], '-depsc', '-painters');
                        saveas(gcf, [TempResultFolder GroupMetaName{iFun} 'CellGroupSLMVsStim'], 'tif');
                
                
                   ResParam.Color=GroupParamNet.FOVMap;
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
                        I0=setdiff(1:NCell,GroupTargetCellAll(GroupTargetCellAll(:,2)==iFun,1));
                        else
                        I0=1:NCell;
                        end
                
                        [Res,r,p]=LuPairRegressPlot_Group(squeeze(rStim(I0,1,1)),statGroupRes(iFun).delta(I0)',Output.NeuroPos3DMeta(I0,4),ResParam,statGroupRes(4).delta(I0));
                        set(gca,'ylim',[-0.1 0.2])
                        if iFun<4
                        title(['Exclude ' GroupLabel{iFun}])
                        end  
                
                        subplotLU(1,4,1,2,P)
                        I0=setdiff(1:NCell,GroupTargetCellAll(:,1));
                        [Res,r,p]=LuPairRegressPlot_Group(squeeze(rStim(I0,1,1)),statGroupRes(iFun).delta(I0)',Output.NeuroPos3DMeta(I0,4),ResParam,statGroupRes(4).delta(I0));
                        set(gca,'YTickLabel',[])
                        set(gca,'ylim',[-0.1 0.2])
                        title(['Exclude 3 groups'])
                
                        subplotLU(1,4,1,3,P)
                        I1=intersect(find(statGroupRes(iFun).delta>0),I0);
                        [Res,r,p]=LuPairRegressPlot_Group(squeeze(rStim(I1,1,1)),statGroupRes(iFun).delta(I1)',Output.NeuroPos3DMeta(I1,4),ResParam,statGroupRes(4).delta(I1));
                        ylabel('')
                        set(gca,'YTickLabel',[])
                                set(gca,'ylim',[-0.1 0.2])
                        title(['Exclude 3 groups'])
                
                         subplotLU(1,4,1,4,P)
                        I1=intersect(find(statGroupRes(iFun).delta<0),I0);
                        [Res,r,p]=LuPairRegressPlot_Group(squeeze(rStim(I1,1,1)),statGroupRes(iFun).delta(I1)',Output.NeuroPos3DMeta(I1,4),ResParam,statGroupRes(4).delta(I1));
                        ylabel('')
                        set(gca,'YTickLabel',[])
                        set(gca,'ylim',[-0.1 0.2])
                        title(['Exclude 3 groups'])
                
                
                
                       clear p

                        papersizePX=[0 0 32 8];
                        set(gcf, 'PaperUnits', 'centimeters');
                        set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
                        LuFontStandard
                        print(gcf, [TempResultFolder GroupMetaName{iFun} 'SLMVsStimRegress.svg'], '-dsvg', '-painters');
                        % print(gcf, [TempResultFolder GroupMetaName{iFun} 'SLMVsStimRegress.eps'], '-depsc', '-painters');
                        saveas(gcf, [TempResultFolder GroupMetaName{iFun} 'SLMVsStimRegress'], 'tif');
                % print(gcf, [TempResultFolder 'SingleSLMVsSpeed.png'], '-dpng', '-painters');
                close all
                end
    
            end
        % end
    end

end





