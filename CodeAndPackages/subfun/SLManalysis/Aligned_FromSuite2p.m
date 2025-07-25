function [AlignedneuralData,SLMInfoTable,CellResponse,statCellRes,TargetResponse,TargetCellResP,TargetCellResR,CellSampleN]=Aligned_FromSuite2p(neuralData,TargetCellList,Suite2pTable,PVpower,PSTHparam)



neuralData=AmpNormalizeRow(neuralData,[0 100]);
SLMInfoTable=Suite2pTable(Suite2pTable.suite2pInd>0,:);
AlignedneuralData=[];

for i=1:size(SLMInfoTable,1)
    I0=SLMInfoTable.suite2pInd(i);
    s1=I0-PSTHparam.PreSLMCal:I0-1;
    s2=I0:I0+PSTHparam.PostSLMCal-1;
    AlignedneuralData(:,:,i)=neuralData(:,[s1 s2])-repmat(nanmean(neuralData(:,s1),2),1,length(s1)+length(s2));
end


CellSampleN=[];
clear TargetResponse CellResponse
TargetCellResP=zeros(length(TargetCellList),length(PVpower))+NaN;
TargetCellResR=zeros(length(TargetCellList),length(PVpower))+NaN;

for iCell=1:length(TargetCellList)
    for iPower=1:length(PVpower)
        I1=find(SLMInfoTable.PointTargetCell==TargetCellList(iCell)&abs(SLMInfoTable.UncagingLaserPower-PVpower(iPower))<0.1);
        CellSampleN(iCell,iPower)=length(I1);

        if length(I1)>1
           % CellResponse{iCell,iPower}=squeeze(AlignedneuralData(:,:,I1));
           % TargetResponse{iCell,iPower}=squeeze(AlignedneuralData(TargetCellList(iCell),:,I1))';

        % elseif length(I1)>1
           CellResponse{iCell,iPower}=nanmean(AlignedneuralData(:,:,I1),3);
           TargetResponse{iCell,iPower}=squeeze(AlignedneuralData(TargetCellList(iCell),:,I1))';
        % else
           preSLMdata=squeeze(AlignedneuralData(:,1:PSTHparam.PreSLMCal,I1));
           postSLMdata=squeeze(AlignedneuralData(:,PSTHparam.PreSLMCal+[1:PSTHparam.TestStepFrame],I1));

           preSLMdata1=squeeze(nanmean(squeeze(AlignedneuralData(:,1:PSTHparam.PreSLMCal,I1)),2));
           postSLMdata1=squeeze(nanmean(squeeze(AlignedneuralData(:,PSTHparam.PreSLMCal+[1:PSTHparam.TestStepFrame],I1)),2));

           preSLMdata2=permute(repmat(preSLMdata1,1,1,PSTHparam.TestStepFrame),[1 3 2]);
           postSLMdata2=postSLMdata-preSLMdata2;

           for jCell=1:size(neuralData,1)
               temp1=preSLMdata(jCell,:,:);
               temp2=postSLMdata(jCell,:,:);
               [~,p(jCell,1),~,t(jCell)]=ttest2(temp2(:),temp1(:));
               temp3(jCell)=mean(temp2(:))-mean(temp1(:));

               % temp1=preSLMdata1(jCell,:);
               % temp2=postSLMdata1(jCell,:);
               % [~,p(jCell,1),~,t(jCell)]=ttest(temp2(:),temp1(:));
               % temp3(jCell)=mean(temp2(:))-mean(temp1(:));
               % 
               % temp2=postSLMdata2(jCell,:,:);
               % [~,p(jCell,1),~,t(jCell)]=ttest(temp2(:),0);
               % temp3(jCell)=mean(temp2(:))-mean(temp1(:));

           end
           statCellRes(iCell,iPower).p=p;
           statCellRes(iCell,iPower).t=t;
           statCellRes(iCell,iPower).delta=temp3;

           TargetCellResP(iCell,iPower)=p(TargetCellList(iCell));
           TargetCellResR(iCell,iPower)=temp3(TargetCellList(iCell));
           clear p t temp3;
        end
        
    end
end
