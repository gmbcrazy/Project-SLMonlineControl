function statGroupRes=OfflineSLM_FOVmeta2NeuroDelta(InputFOVData, Par,PSTHparam,SaveP1)

StructToVars(Par);
StructToVars(InputFOVData);

GroupMetaName = [GroupLabel {'FakeSLM'}];
GroupMetaColor = [GroupColor; PowerZeroColor];
NCell = size(InputFOVData.NeuroPos3DMeta,1);
FunNum = length(GroupList) + 1;
TimBinFrame = -PSTHparam.PreSLMCal:PSTHparam.PostSLMCal-1;

TestStepFrame=PSTHparam.TestStepFrame;
        % clear statGroupRes

for iVol = 1:length(VolOut)
    for iAwake = 1:length(AwakeState)

        % clear GroupResponse GroupSpeed statGroupRes GroupSampleN

        statGroupRes(iVol,iAwake).p         = ones(NCell, FunNum);
        statGroupRes(iVol,iAwake).delta     = zeros(NCell, FunNum);
        statGroupRes(iVol,iAwake).speeddelta= zeros(NCell, FunNum);
        statGroupRes(iVol,iAwake).GroupSampleN = zeros(NCell, FunNum);
        statGroupRes(iVol,iAwake).GroupSpeed = nan(NCell, length(TimBinFrame), FunNum);
        statGroupRes(iVol,iAwake).GroupResponse = nan(NCell, length(TimBinFrame), FunNum);
        statGroupRes(iVol,iAwake).GroupSpeedTrial = {};


       for iFun = 1:length(GroupList)
           I1All=find(AlignedInfoTable.Group==GroupList(iFun)&AlignedInfoTable.PowerZero==0&AlignedInfoTable.VolOut==VolOut(iVol)&AlignedInfoTable.AwakeState==AwakeState(iAwake));
            % GroupSampleN(iFun)=length(I1);

           preSLMSpeed=mean(squeeze(AlignedSpeedMeta(1:PSTHparam.PreSLMCal,I1All)),1);
           postSLMSpeed=mean(squeeze(AlignedSpeedMeta(PSTHparam.PreSLMCal+[1:TestStepFrame],I1All)),1);

           % % if iAwake==1
           % % I1All=I1All(find(abs(postSLMSpeed-preSLMSpeed)<=PostPreDiffSpeedTh));
           % % end
           % % 
           % % GroupSpeed{iFun}=AlignedSpeedMeta(:,I1All)';
           % % GroupSampleN(iFun)=length(I1All);
           % iFun
           SpeedGroup{iFun}=[];
           for iFOV = 1:length(InputFOVData.AlignedNData)
               % iFOV
                AlignedInfoTableFOV=InputFOVData.AlignedInfoTableFOV{iFOV};
                I1=find(AlignedInfoTableFOV.Group==GroupList(iFun)&AlignedInfoTableFOV.PowerZero==0&AlignedInfoTableFOV.VolOut==VolOut(iVol)&AlignedInfoTableFOV.AwakeState==AwakeState(iAwake));
                AlignedtempNData = InputFOVData.AlignedNData{iFOV};
                iFOVneuro=find(InputFOVData.NeuroPos3DMeta(:,4)==iFOV);
                iFOVEvent = find(AlignedInfoTable.iFOV==iFOV);

                AlignedSpeedFOV=AlignedSpeedMeta(:,iFOVEvent);
                preSLMSpeed=mean(squeeze(AlignedSpeedFOV(1:PSTHparam.PreSLMCal,I1)),1);
                postSLMSpeed=mean(squeeze(AlignedSpeedFOV(PSTHparam.PreSLMCal+[1:TestStepFrame],I1)),1);


                if iAwake==1
                   I1temp=find(abs(postSLMSpeed-preSLMSpeed)<=PostPreDiffSpeedTh);
                   I1=I1(I1temp);
                end
                preSLMSpeed=preSLMSpeed(I1temp);
                postSLMSpeed=postSLMSpeed(I1temp);
                SpeedGroup{iFun}=[SpeedGroup{iFun};[preSLMSpeed(:) postSLMSpeed(:)]];




                [~,pSpeed(iFun, iFOV),~,statSpeed(iFun, iFOV)]=ttest(postSLMSpeed,preSLMSpeed);
                SpeedDelta = mean(postSLMSpeed) - mean(preSLMSpeed);
                % size(statGroupRes(iVol,iAwake).speeddelta)
                statGroupRes(iVol,iAwake).speeddelta(iFOVneuro,iFun)=SpeedDelta;

                if length(I1)==1
                   % statGroupRes(iVol,iAwake).GroupResponse(iFOVneuro,:,iFun)=squeeze(AlignedtempNData(:,:,I1));
                   statGroupRes(iVol,iAwake).GroupSpeed(iFOVneuro,:,iFun)=repmat(AlignedSpeedMeta(:,I1)',length(iFOVneuro), 1);
                   statGroupRes(iVol,iAwake).GroupResponse(iFOVneuro,:,iFun)=squeeze(AlignedtempNData(:,:,I1));
                   statGroupRes(iVol,iAwake).GroupSampleN(iFOVneuro,iFun)=length(I1);

                   preSLMdata=squeeze(AlignedtempNData(:,1:PSTHparam.PreSLMCal,I1));
                   postSLMdata=squeeze(AlignedtempNData(:,PSTHparam.PreSLMCal+[1:TestStepFrame],I1));
                   for jCell=1:sum(InputFOVData.NeuroPos3DMeta(:,4)==iFOV)
                       temp1=preSLMdata(jCell,:,:);
                       temp2=postSLMdata(jCell,:,:);
                       [~,p(jCell,1),~,t(jCell)]=ttest(temp2(:),temp1(:));
                       temp3(jCell)=mean(temp2(:))-mean(temp1(:));
                   end
                elseif length(I1)>1
                   % GroupSpeed{iFun}=AlignedSpeedMeta(:,I1)';
                   statGroupRes(iVol,iAwake).GroupSpeed(iFOVneuro,:,iFun)=repmat(nanmean(AlignedSpeedMeta(:,I1),2)',length(iFOVneuro), 1);
                   statGroupRes(iVol,iAwake).GroupResponse(iFOVneuro,:,iFun)=nanmean(AlignedtempNData(:,:,I1),3);  
                   preSLMdata=squeeze(AlignedtempNData(:,1:PSTHparam.PreSLMCal,I1));
                   postSLMdata=squeeze(AlignedtempNData(:,PSTHparam.PreSLMCal+[1:TestStepFrame],I1));
                   for jCell=1:sum(InputFOVData.NeuroPos3DMeta(:,4)==iFOV)
                       temp1=squeeze(nanmean(preSLMdata(jCell,:,:),2));
                       temp2=squeeze(nanmean(postSLMdata(jCell,:,:),2));
                       [~,p(jCell,1),~,t(jCell)]=ttest(temp2(:),temp1(:));
                       temp3(jCell)=mean(temp2(:))-mean(temp1(:));
                   end
                   statGroupRes(iVol,iAwake).p(iFOVneuro,iFun)=p;
                   % statGroupRes(iFun).t=t;
                   statGroupRes(iVol,iAwake).delta(iFOVneuro,iFun)=temp3;
                   statGroupRes(iVol,iAwake).GroupSampleN(iFOVneuro,iFun)=length(I1);
                   clear p t temp3;
        
        
                else
                   % statGroupRes(iFun).p=zeros(size(iscell))+1;
                   % statGroupRes(iFun).delta=zeros(size(iscell));

                end

        end


            
        end
    
    

        % statGroupRes(iFun).p=zeros(NCell,1)+1;
        % statGroupRes(iFun).delta=zeros(NCell,1);
        % GroupResponse{iFun}=zeros(NCell,length(TimBinFrame));

        iFun=iFun+1;
%%Add Zero power as addtional group
        I1All=find(AlignedInfoTable.PowerZero==1&AlignedInfoTable.VolOut==VolOut(iVol)&AlignedInfoTable.AwakeState==AwakeState(iAwake));
        preSLMSpeed=mean(squeeze(AlignedSpeedMeta(1:PSTHparam.PreSLMCal,I1All)),1);
        postSLMSpeed=mean(squeeze(AlignedSpeedMeta(PSTHparam.PreSLMCal+[1:TestStepFrame],I1All)),1);
        % if iAwake==1
        %    I1All=I1All(find(abs(postSLMSpeed-preSLMSpeed)<=PostPreDiffSpeedTh));
        % end
        % GroupSpeed{iFun}=AlignedSpeedMeta(:,I1All)';

        % GroupSampleN(iFun)=length(I1All);
        SpeedGroup{iFun}=[];
        for iFOV = 1:length(InputFOVData.AlignedNData)
                AlignedInfoTableFOV=InputFOVData.AlignedInfoTableFOV{iFOV};
                I1=find(AlignedInfoTableFOV.PowerZero==1&AlignedInfoTableFOV.VolOut==VolOut(iVol)&AlignedInfoTableFOV.AwakeState==AwakeState(iAwake));
                AlignedtempNData = InputFOVData.AlignedNData{iFOV};
                iFOVneuro=find(InputFOVData.NeuroPos3DMeta(:,4)==iFOV);


                iFOVEvent = find(AlignedInfoTable.iFOV==iFOV);
                AlignedSpeedFOV=AlignedSpeedMeta(:,iFOVEvent);
                preSLMSpeed=mean(squeeze(AlignedSpeedFOV(1:PSTHparam.PreSLMCal,I1)),1);
                postSLMSpeed=mean(squeeze(AlignedSpeedFOV(PSTHparam.PreSLMCal+[1:TestStepFrame],I1)),1);



                if iAwake==1
                    I1temp=find(abs(postSLMSpeed-preSLMSpeed)<=PostPreDiffSpeedTh);
                    I1=I1(I1temp);
                end
                preSLMSpeed=preSLMSpeed(I1temp);
                postSLMSpeed=postSLMSpeed(I1temp);
                SpeedGroup{iFun}=[SpeedGroup{iFun};[preSLMSpeed(:) postSLMSpeed(:)]];



                [~,pSpeed(iFun,iFOV),~,statSpeed(iFun,iFOV)]=ttest(postSLMSpeed,preSLMSpeed);
                SpeedDelta = mean(postSLMSpeed) - mean(preSLMSpeed);
                statGroupRes(iVol,iAwake).speeddelta(iFOVneuro,iFun)=SpeedDelta;






            if length(I1)==1
               % GroupResponse{iFun}=squeeze(AlignedtempNData(:,:,I1));
               % GroupSpeed{iFun}=AlignedSpeedMeta(:,I1)';
               statGroupRes(iVol,iAwake).GroupSpeed(iFOVneuro,:,iFun)=repmat(AlignedSpeedMeta(:,I1)',length(iFOVneuro),1);
               statGroupRes(iVol,iAwake).GroupResponse(iFOVneuro,:,iFun)=squeeze(AlignedtempNData(:,:,I1));
               statGroupRes(iVol,iAwake).GroupSampleN(iFOVneuro,iFun)=length(I1);

            elseif length(I1)>1
               % GroupResponse{iFun}=nanmean(AlignedtempNData(:,:,I1),3);    
               statGroupRes(iVol,iAwake).GroupResponse(iFOVneuro,:,iFun)=nanmean(AlignedtempNData(:,:,I1),3);  
               % GroupSpeed{iFun}=AlignedSpeedMeta(:,I1)';
               % statGroupRes(iVol,iAwake).GroupSpeed(iFOVneuro,:,iFun)=nanmean(AlignedSpeedMeta(:,I1),2);
               statGroupRes(iVol,iAwake).GroupSpeed(iFOVneuro,:,iFun)=repmat(nanmean(AlignedSpeedMeta(:,I1),2)',length(iFOVneuro), 1);

               preSLMdata=squeeze(AlignedtempNData(:,1:PSTHparam.PreSLMCal,I1));
               postSLMdata=squeeze(AlignedtempNData(:,PSTHparam.PreSLMCal+[1:TestStepFrame],I1));
               for jCell=1:sum(InputFOVData.NeuroPos3DMeta(:,4)==iFOV)

                    temp1=squeeze(nanmean(preSLMdata(jCell,:,:),2));
                    temp2=squeeze(nanmean(postSLMdata(jCell,:,:),2));
                    [~,p(jCell,1),~,t(jCell)]=ttest(temp2(:),temp1(:));
                    temp3(jCell)=mean(temp2(:))-mean(temp1(:));

               end
               statGroupRes(iVol,iAwake).p(iFOVneuro,iFun)=p;
               % statGroupRes(iFun).t=t;
               statGroupRes(iVol,iAwake).delta(iFOVneuro,iFun)=temp3;
               statGroupRes(iVol,iAwake).GroupSampleN(iFOVneuro,iFun)=length(I1);

               clear p t temp3;
    
    
            else
               % statGroupRes(iFun).p=zeros(size(iscell))+1;
               % statGroupRes(iFun).delta=zeros(size(iscell));

            end

        end
 %%Add Zero power as addtional group
    
        statGroupRes(iVol,iAwake).deltaNorm=statGroupRes(iVol,iAwake).delta(:,length(GroupList))-repmat(statGroupRes(iVol,iAwake).delta(:,length(GroupList)+1),1,length(GroupList));

        statGroupRes(iVol,iAwake).GroupSpeedTrial=SpeedGroup;
        clear SpeedGroup

    end
end



