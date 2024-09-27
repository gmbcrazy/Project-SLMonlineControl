function [UpdateXml, SLMTable, PointsTest, XMLparam, InfoListByLaser]=step3Fun_NextSLMtest(SLMRes,sampleN,ROIparam,XMLparam,SLMTestParam,SLMTable)

TerminalTrialN=SLMTestParam.TerminalTrialN;  %Trials # to define SLM responsive cells
ExcludeTrialN=SLMTestParam.ExcludeTrialN;   %Trials # to define SLM non-responsive cells


PointsTest=[];
UpdateXml=0;

% Potential next round of Non-targets were chosen.
NextRoundID=randperm(XMLparam.TotalRounds,1);


for iLaser=1:length(ROIparam.LaserPower)
    % [~,iLaser]=ismember(ROIparam.LaserPower(iLaser),confSet.UncagingLaserPower);
    if iLaser==1
       CheckedPoints = [];
    else
       CheckedPoints=find(sum(find(SLMRes(:,1:iLaser-1)==1&sampleN(:,1:iLaser-1)>=TerminalTrialN),2)>0);
    end
    

    SLMList{iLaser}=find(SLMRes(:,iLaser)==1&sampleN(:,iLaser)>=TerminalTrialN);
    PriorityConfirmList{iLaser}=find(SLMRes(:,iLaser)==1&sampleN(:,iLaser)<TerminalTrialN);
    ExcludeList{iLaser}=find(SLMRes(:,iLaser)==0&sampleN(:,iLaser)>=ExcludeTrialN);
    UnkownConfirmList{iLaser}=find(SLMRes(:,iLaser)==0&sampleN(:,iLaser)<ExcludeTrialN);

    SLMList{iLaser}=setdiff(SLMList{iLaser},CheckedPoints);                            %Exclude previous responsive cells for current level 
    PriorityConfirmList{iLaser}=setdiff(PriorityConfirmList{iLaser},CheckedPoints);    %Exclude previous responsive cells for current level 
    ExcludeList{iLaser}=setdiff(ExcludeList{iLaser},CheckedPoints);                    %Exclude previous responsive cells for current level
    UnkownConfirmList{iLaser}=setdiff(UnkownConfirmList{iLaser},CheckedPoints);        %Exclude previous responsive cells for current level 




    InfoListByLaser(iLaser).SLMList=SLMList{iLaser};
    InfoListByLaser(iLaser).PriorityConfirmList=PriorityConfirmList{iLaser};
    InfoListByLaser(iLaser).ExcludeList=ExcludeList{iLaser};
    InfoListByLaser(iLaser).UnkownConfirmList=UnkownConfirmList{iLaser};
    InfoListByLaser(iLaser).LaserPower=ROIparam.LaserPower(iLaser);


    SLMtempTable=SLMTable;
    SLMtempTable(:,3)=NaN;
    SLMtempTable(SLMList{iLaser},3)=ROIparam.LaserPower(iLaser);

    SLMTable(:,2)=min(SLMtempTable(:,[2 3]),[],2);
   %%
   % if (length(PriorityConfirmList{iLaser})+PriorityConfirmList{iLaser})/
   % 
   % 
   % end

   %% 
   if ~isempty(SLMList{iLaser})
       disp([num2str(length(SLMList{iLaser})) ' Priority Points respond to laser power' num2str(ROIparam.LaserPower(iLaser)) ': ' num2str(SLMList{iLaser}') ]);
   end
   %%
   if ~isempty(ExcludeList{iLaser})
      disp([num2str(length(ExcludeList{iLaser})) ' Points were excluded for laser power' num2str(ROIparam.LaserPower(iLaser)),' requiring higher power: ' num2str(ExcludeList{iLaser}') ]);
   end

    %% First test the MP point with response but trial # is not enough
    if ~isempty(PriorityConfirmList{iLaser})&&UpdateXml==0
       UpdateXml=1;
       PointsTest=PriorityConfirmList{iLaser};
       XMLparam.Laser=ROIparam.LaserPower(iLaser);
       XMLparam.RoundID=NextRoundID;
       disp([num2str(length(PriorityConfirmList{iLaser})) ' Priority Points need more test with laser power' num2str(ROIparam.LaserPower(iLaser)) ':' num2str(PriorityConfirmList{iLaser}') ]);
       % continue
       % break
    end

    %% Then, test the MP point without response while trial # is not enough
    if ~isempty(PriorityConfirmList{iLaser})&&UpdateXml==0
       UpdateXml=1;
       PointsTest= UnkownConfirmList{iLaser};
       XMLparam.RoundID=NextRoundID;
       XMLparam.Laser=ROIparam.LaserPower(iLaser);
       disp([num2str(length(UnkownConfirmList{iLaser})) ' Unkown Points need more test with laser power' num2str(ROIparam.LaserPower(iLaser)) ': ' num2str(UnkownConfirmList{iLaser}') ]);
       % break
    end

    if isempty(PointsTest)&&(~isempty(ExcludeList{iLaser}))
       I1=find(SLMTestParam.AllLaserPower>ROIparam.LaserPower(iLaser));
       if isempty(I1)&&UpdateXml==0
          PointsTest= UnkownConfirmList{iLaser};      
          XMLparam.Laser=SLMTestParam.AllLaserPower(min(I1));
          XMLparam.RoundID=NextRoundID;
          UpdateXml=1;
          disp('Warning: power is updated to higher than ROIparam levels');
       else
          disp(['Warning: No higher power levels could be test with current Unconfirmed Points: ' num2str(PointsTest')]);
       end
     
    end

end

if UpdateXml==1
   disp(['Next Points will be test: ' num2str(PointsTest') ' at power level ' num2str(XMLparam.Laser)]);
end




