function [PointLaserPair,ResPointLaser] = SelectPointsForTesting_v4(SLMRes, sampleN, SLMTestParam, ZTestNum)


%%SelectPointsForTesting_v4 is modified from v3, Lu Zhang, 01282025
%where v3 extract next testing PointList with possible different laser level.
%where v4 extract testing PointList all under the same laser level


    % Get the size of SLMRes and sampleN
    [N, L] = size(SLMRes);    %%N is the total PossibleTest Points, L is number of laser level for testing
    
    % N =10;
    % L=3;
    % Generate the PointsListAll and LaserListAll in a vectorized manner
    [iPoints, iTrials, iLasers] = ndgrid(1:N, 1:SLMTestParam.TerminalTrialN, 1:L);


    PointsListAll = iPoints(:);
    LaserListAll = iLasers(:);
    TrialListAll = iTrials(:);

    PointsListAllRaw=PointsListAll;
    LaserListAllRaw=LaserListAll;
    TrialListAllRaw=TrialListAll;

    
    ResPoint=find(sum(SLMRes==1&sampleN>=SLMTestParam.TerminalTrialN,2)>=1); %responsive cell with enough trials tested, no need to test further more
    DelI=find(ismember(PointsListAll,ResPoint)>=1);


    [ResPointID,ResPointLaser]=find(SLMRes==1&sampleN>=SLMTestParam.TerminalTrialN); %responsive cell with enough trials tested, no need to test further more
    [ResPointID, firstIdx] = unique(ResPointID, 'stable');
    ResPointLaser=ResPointLaser(firstIdx);
    ResPointLaser=[ResPointID ResPointLaser];
    % ResPoint
    % PointsListAll(Del0)=[];
    % LaserListAll(Del0)=[];
    % TrialListAll(Del0)=[];



    for iLaser=1:size(sampleN,2)
        NonResPoint=find(SLMRes(:,iLaser)==0&sampleN(:,iLaser)>=SLMTestParam.ExcludeTrialN);  %%No responsive cell with enough trials tested
        DelI1=ismember(PointsListAll,NonResPoint)&LaserListAll<=iLaser;                       %%No need to test current laser and lower laser

        NonResPoint=find(SLMRes(:,iLaser)==0&sampleN(:,iLaser)<SLMTestParam.ExcludeTrialN);   %%No responsive cell needs be further tested
        DelI2=ismember(PointsListAll,NonResPoint)&LaserListAll>iLaser;                        %%Test current laser, not the higher         

        % NonResPoint=find(SLMRes(:,iLaser)==0&sampleN(:,iLaser)<SLMTestParam.ExcludeTrialN);  %%No responsive cell needs be further tested        
        % DelI3=ismember(PointsListAll,NonResPoint)&LaserListAll==iLaser&TrialListAll>SLMTestParam.ExcludeTrialN;    %%Test current laser, but not too many trials within a Tseries         
%         DelI3=ismember(PointsListAll,NonResPoint)&LaserListAll==iLaser&TrialListAll>=SLMTestParam.TerminalTrialN;    %%Test current laser, but not too many trials within a Tseries         
   
        ResPoint=find(SLMRes(:,iLaser)==1&sampleN(:,iLaser)<SLMTestParam.TerminalTrialN);  %%Responsive cell needs be further tested        
        DelI4=[];
        for iP=1:length(ResPoint)
            DelItemp1=PointsListAll==ResPoint(iP)&LaserListAll>iLaser;    %%Test current laser
%             DelItemp2=PointsListAll==ResPoint(iP)&LaserListAll==iLaser&TrialListAll>SLMTestParam.TerminalTrialN-sampleN(ResPoint(iP),iLaser)+2;    %%but not too many trials within a Tseries                 
%             DelI4=union(DelI4,find(DelItemp1+DelItemp2>=1));
            DelI4=union(DelI4,find(DelItemp1>=1));
          
        end

        DelTemp1=union(find(DelI1+DelI2>=1),DelI4);
        DelI=union(DelI,DelTemp1);
        
    end
 
    PointsListAll(DelI)=[];
    LaserListAll(DelI)=[];
    TrialListAll(DelI)=[];

    DelI=[];
    for iLaser=1:size(sampleN,2)
        for iPoint=1:N
            tempI=find(PointsListAll==iPoint&TrialListAll<=sampleN(iPoint,iLaser)&LaserListAll==iLaser);
            DelI=union(DelI,tempI(:));
        end
    end

    PointsListAll(DelI)=[];
    LaserListAll(DelI)=[];
    TrialListAll(DelI)=[];


    if isempty(PointsListAll)
       PointLaserPair=[];
       return
    end

    %%Only test minal laser level needs to tested
    minLaser=min(LaserListAll);
    DelI=find(LaserListAll>minLaser);
    PointsListAll(DelI)=[];
    LaserListAll(DelI)=[];
    TrialListAll(DelI)=[];



%     [TrialListAll,sortI]=sort(TrialListAll);
%     PointsListAll=PointsListAll(sortI);
%     LaserListAll=LaserListAll(sortI);


    %%There are still many points and power levels to test;
    if ZTestNum<=length(PointsListAll)
       TestPoints=PointsListAll(1:ZTestNum);
       TestLaserLevels=LaserListAll(1:ZTestNum);
       TestTrial=TrialListAll(1:ZTestNum);
       PointLaserPair=[TestPoints TestLaserLevels];
       return
    end
 
    %%There are more than enough points to test, add some points;
    PointLaserPair=unique([PointsListAll LaserListAll],'rows');

    I1=find(ResPointLaser(:,2)==minLaser); %Points with current laser works.
    ResPointLaser=ResPointLaser(I1,:);

    if length(ResPointLaser)>=3      
       AddPair=ResPointLaser;  %%existing point works, add them simply increase their samples
    else
       AddP=setdiff([1:N]',PointLaserPair(:,1));  %%almost no existing point works, add them with current lowest power level which do not work anyway
       AddPair=[AddP ones(length(AddP),1)*minLaser];
    end
    N0=size(PointLaserPair,1);
    N1=size(AddPair,1);


       FinalPairtemp=[];
       if N0<=2
          while size(FinalPairtemp,1)<ZTestNum
                 addI1=randperm(N1,min(N1,2));    %%add at most 2 points from addpair to ensure most test are done for non tested points and power pair.
                 FinalPairtemp=[FinalPairtemp;PointLaserPair;AddPair(addI1,:)];
          end
       else
          while size(FinalPairtemp,1)<ZTestNum
                 % addI1=randperm(N1,min(N1,2));    
                 FinalPairtemp=[FinalPairtemp;PointLaserPair]; %%simply repeated test the exising points and laser level needs to be tested.
          end
 
       end
       PointLaserPair=FinalPairtemp(1:ZTestNum,:);





    % figure;hold on;
    % plot(PointLaserPair(:,1),'r.-');
    % plot(PointLaserPair(:,2),'g^');
    % % plot(TestTrial,'bo');
    % legend({'Point','laser'})
    % 

end

    


