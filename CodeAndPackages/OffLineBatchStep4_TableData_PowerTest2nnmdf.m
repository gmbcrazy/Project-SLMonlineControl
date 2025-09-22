
load('\\nimhlabstore1.nimh.nih.gov\UFNC\FNC2\Zhang\Projects\Project-LocalProcessing\Step4\awakeRefSpon\GroupSLM16Sessions\FunCon\FC_pth05\deltaFVecTable.mat');

load('\\nimhlabstore1.nimh.nih.gov\UFNC\FNC2\Zhang\Projects\Project-LocalProcessing\Step4\awakeRefSpon\GroupSLM16Sessions\NonNegMatFac\FOV\deltaF29-Aug-2025nnmfResult.mat')


PowerTestTBL=groupsummaryBack2OldNames(tblFOVclean,tblFOVcleanavg,'mean');

[G, sessionVals, TargetcellVals] = findgroups(PowerTestTBL.Session, PowerTestTBL.PointTargetCell);

UniqueG=unique(G);
UniqueS=unique(sessionVals);
tblBridge=table(sessionVals,TargetcellVals,'VariableNames',{'Session','PointTargetCell'});
newtbl5 = groupsummary(PowerTestTBL, {'Session', 'Cell'}, 'mean');

newtblScore = groupsummary(PowerTestTBL, {'Session', 'PointTargetCell'}, 'mean');
PowerScoreTBL=groupsummaryBack2OldNames(PowerTestTBL,newtblScore,'mean');


PowerScore=[];

for iG=1:length(UniqueG)
    iFOV=find(UniqueS==sessionVals(iG));
    jFOV=UniqueS(iFOV);
    % TargetcellVals(iG)
    % if ismember(iFOV,[5,6])
    %    continue;
    % end
    TargetCellTemp=TargetcellVals(iG);
    % A1=find(newtbl5.Session==sessionVals(iG));
    % B1=find(newtbl5.mean_ActivateI(A1)>CoActITh);
    % B2=find(newtbl5.mean_ActivateIPos(A1)>CoActITh);
    % B3=find(newtbl5.mean_ActivateINeg(A1)>CoActITh);
    % testtbl=PowerTestTBL(G==UniqueG(iG),:);

    AllInd=find(G==UniqueG(iG));
    AllIndSub=find(G==UniqueG(iG)&PowerTestTBL.TargetCell==0);

    tempVec=PowerTestTBL.Response(AllInd);
    tempVecSub=PowerTestTBL.Response(AllIndSub);

    % % tempActI{iFOV}{end+1} = PowerTestTBL.ActivateI(G==UniqueG(iG));
    % % tempActIPos{iFOV}{end+1} = PowerTestTBL.ActivateIPos(G==UniqueG(iG));
    % % tempActINeg{iFOV}{end+1} = PowerTestTBL.ActivateINeg(G==UniqueG(iG));
    % % 
    % % 
    % % tempVecAct{iFOV}(:,end+1)=tempVec{iFOV}(B1,end);
    % % tempVecActP{iFOV}(:,end+1)=tempVec{iFOV}(B2,end);
    % % tempVecActN{iFOV}(:,end+1)=tempVec{iFOV}(B3,end);
    % PowerScore=lsqnonneg(FOVres(iFOV).W(:,FOVres(iFOV).ComI), tempVec);
    % PowerScoreNontarget=lsqnonneg(FOVres(iFOV).WNontarget(:,FOVres(iFOV).ComI), tempVecSub);

    pinvW = pinv(FOVres(jFOV).W(:,FOVres(jFOV).ComI));
    pinvW_sub = pinv(FOVres(jFOV).W_sub(:,FOVres(jFOV).ComI_sub));


    Score1 = (pinvW * tempVec)'; 
    Score_sub1 = (pinvW_sub * tempVecSub)';

    PowerScore=[PowerScore;Score1 Score_sub1 jFOV TargetCellTemp];
end

PowerScore=array2table(PowerScore,'VariableNames',{'SpeedScore','StimScore','SpeedScoreNontarget','StimScoreNontarget','Session','PointTargetCell'});


sumT1 = join(PowerScore, PowerScoreTBL, 'Keys', {'Session','PointTargetCell'});
sumT1(:,{'TargetSpeedR','TargetStimR'})=[];


FactorVariGroupComp={'SpeedScore','SpeedScoreNontarget','StimScore','StimScoreNontarget'};
figure;
for jFF=1:length(FactorVariGroupComp)
    subplot(1,length(FactorVariGroupComp),jFF)
    clear tempdata
    for igroup=1:3
        tempdata{igroup}=table2array(sumT1(sumT1.PointTargetCellGroup==igroup,FactorVariGroupComp{jFF}));
    %     histPlotLU(sumT1.VecR(sumT1.PointTargetCellGroup==igroup),DisX,ProcessPar.GroupColor(igroup,:),0.5);
    % hold on;
    end
   GroupPair.CorrName='fdr';
   GroupPair.Q=0.1;
   GroupPair.Pair=[];
   GroupPair.SignY=0.02;
   GroupPair.Plot=1;
   GroupPair.Std=0;      %%%%%%%%%using standard deviation as errorbar
   GroupPair.SamplePlot=1; %%%%%%%%%Plot Individual Sample Point
   GroupPair.SamplePairedPlot=0; %%%%%%%%%Dash line for paired comparison sample
   GroupPair.LimY=[-0.003 GroupPair.SignY*1.2];
   GroupPair.Marker={'o'};   % Datatype=varargin{1};  %%%%%%Datatype=0, non-paired data; 1 for paired data
   GroupPair.GroupName={'L','S','N'}
   GroupPairViolin=GroupPair;
   GroupPairViolin.Test='Ranktest';
   GroupPairViolin.ViolinLR=[0 0 0]
   % subplot(1,2,2)
   % stats=ErrorBarPlotLU(1:3,tempdata,[],ProcessPar.GroupColor,2,0,'VecTrialStabilityGroup.txt',GroupPair,[1 2 3]);
   % stats=ErrorViolinHalf(1:3,tempdata,ProcessPar.GroupColor,0,[SaveFunConPth  NDataName{iData} FactorVariGroupComp{jFF} 'Group.txt'],GroupPairViolin,[1 2 3])
   stats(jFF)=ErrorViolinHalf(1:3,tempdata,ProcessPar.GroupColor,0,[],GroupPairViolin,[1 2 3])

   xticks(1:3)
   xticklabels(GroupPairViolin.GroupName)
   ylabel(FactorVariGroupComp{jFF})
   set(gca,'ylim',[-0.00001 0.00002]);

end



Factor={'PointSpeedR','SpeedR','PointStimR','StimR','VecR','ActVecR'};

close all
ParRegress.Color=ProcessPar.GroupColor;
ParRegress.Marker='s';
ParRegress.MarkerSize=10;
ParRegress.Rtype='Spearman';
ParRegress.xLim=[0 0.1];
ParRegress.yLim=[0 0.2];
% ParRegress.xLabel='StimR';
% ParRegress.yLabel=Predictor{iPP};
% [OutPut,r,p]=LuPairRegressPlot_Group(newtbl3.StimR,newtbl3.ActivateI,newtbl3.PointTargetCellGroup,ParRegress)    
ParRegress.xLim=[];
ParRegress.yLim=[];


SubplotParam=[0.05 0.03 0.1 0.2 0.03 0.1]

for iFF=1:length(Factor)
    figure;
    for iPP=1:length(FactorVariGroupComp)
        ParRegress.yLabel=FactorVariGroupComp{iPP};
        ParRegress.xLabel=Factor{iFF};
        ParRegress.xLim=[];
        ParRegress.yLim=[];

        data1=table2array(sumT1(:,Factor{iFF}));
        data2=table2array(sumT1(:,FactorVariGroupComp{iPP}));

        dataCov=table2array(sumT1(:,setdiff(Factor,Factor{iFF})));

        ValidI=(~isnan(data1))|(~isnan(data2));
        data1=data1(ValidI);
        data2=data2(ValidI);
        dataCov=dataCov(ValidI,:);

        h(1)=subplotLU(2,length(FactorVariGroupComp),1,iPP,SubplotParam);
        [OutPut,r,p]=LuPairRegressPlot_Group(data1,data2,sumT1.PointTargetCellGroup(ValidI),ParRegress)
        xlabel('');
        % a=title(FactorVariGroupComp{iPP})

        h(1)=subplotLU(2,length(FactorVariGroupComp),2,iPP,SubplotParam);
        [B,BINT,R,RINT,STATS]=regress(data2,[ones(length(dataCov),1) dataCov]);
        ParRegress.yLim=[];
        [OutPut,r,p]=LuPairRegressPlot_Group(data1,R,sumT1.PointTargetCellGroup(ValidI),ParRegress);

            
    end

    papersizePX=[0 0 length(FactorVariGroupComp)*4.5 10];
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
    % print(gcf, [SaveFunConPth  NDataName{iData} Factor{iFF} 'MixedGroup.svg'], '-dsvg', '-painters');
    % print(gcf, [SaveFunConPth  NDataName{iData} Factor{iFF} 'MixedGroup.tif'], '-dtiffn', '-painters');
    % 

end

SubplotParam=[0.05 0.03 0.025 0.05 0.025 0.025]
ParRegressSub=ParRegress;
for iFF=1:length(Factor)
    figure;
    for iPP=1:length(FactorVariGroupComp)
        ParRegressSub.yLabel=FactorVariGroupComp{iPP};
        ParRegressSub.xLabel=Factor{iFF};

        for igroup=1:3
        ParRegressSub.Color=ParRegress.Color(igroup,:);
        ParRegressSub.xLim=[];
        ParRegressSub.yLim=[];
        subI=find(sumT1.PointTargetCellGroup==igroup);
        data1=table2array(sumT1(subI,Factor{iFF}));
        data2=table2array(sumT1(subI,FactorVariGroupComp{iPP}));

        dataCov=table2array(sumT1(subI,setdiff(Factor,Factor{iFF})));

        ValidI=(~isnan(data1))|(~isnan(data2));
        data1=data1(ValidI);
        data2=data2(ValidI);
        dataCov=dataCov(ValidI,:);


        h(1)=subplotLU(6,length(FactorVariGroupComp),igroup,iPP,SubplotParam);
        [OutPut,r,p]=LuPairRegressPlot_Group(data1,data2,ones(size(data1)),ParRegressSub);

        xlabel('');

        % a=title(FactorVariGroupComp{iPP})

        h(1)=subplotLU(6,length(FactorVariGroupComp),igroup+3,iPP,SubplotParam);
        [B,BINT,R,RINT,STATS]=regress(data2,[ones(length(dataCov),1) dataCov]);
        ParRegressSub.yLim=[min(R) max(R)];
        [OutPut,r,p]=LuPairRegressPlot_Group(data1,R,ones(size(data1)),ParRegressSub);
        if igroup~=3
        xlabel('');
        end
        end
    end

    papersizePX=[0 0 length(FactorVariGroupComp)*4.5 5*6];
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
    % print(gcf, [SaveFunConPth  NDataName{iData} Factor{iFF} 'SepGroup.svg'], '-dsvg', '-painters');
    % print(gcf, [SaveFunConPth  NDataName{iData} Factor{iFF} 'SepdGroup.tif'], '-dtiffn', '-painters');


end
close all
