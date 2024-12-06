function PointLaserPair = SelectPointsForTesting_v3(SLMRes, sampleN, SLMTestParam, ZTestNum)
    % Get the size of SLMRes and sampleN
    [N, L] = size(SLMRes);    %%N is the total PossibleTest Points, L is number of laser level for testing
    
    % N =10;
    % L=3;
    % Generate the PointsListAll and LaserListAll in a vectorized manner
    [iPoints, iTrials, iLasers] = ndgrid(1:N, 1:SLMTestParam.TerminalTrialN, 1:L);
    PointsListAll = iPoints(:);
    LaserListAll = iLasers(:);
    TrialListAll = iTrials(:);

    
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
        DelI1=ismember(PointsListAll,NonResPoint)&LaserListAll==iLaser;
        % PointsListAll(DelI1)=[];
        % LaserListAll(DelI1)=[];
        % TrialListAll(DelI1)=[];

        NonResPoint=find(SLMRes(:,iLaser)==0&sampleN(:,iLaser)<SLMTestParam.ExcludeTrialN);  %%No responsive cell needs be further tested
        DelI2=ismember(PointsListAll,NonResPoint)&LaserListAll>iLaser;                        %%Test current laser, not the higher         
        % PointsListAll(DelI2)=[];
        % LaserListAll(DelI2)=[];
        % TrialListAll(DelI2)=[];

        NonResPoint=find(SLMRes(:,iLaser)==0&sampleN(:,iLaser)<SLMTestParam.ExcludeTrialN);  %%No responsive cell needs be further tested        
        % DelI3=ismember(PointsListAll,NonResPoint)&LaserListAll==iLaser&TrialListAll>SLMTestParam.ExcludeTrialN;    %%Test current laser, but not too many trials within a Tseries         
        DelI3=ismember(PointsListAll,NonResPoint)&LaserListAll==iLaser&TrialListAll>=SLMTestParam.TerminalTrialN;    %%Test current laser, but not too many trials within a Tseries         

        % PointsListAll(DelI3)=[];
        % LaserListAll(DelI3)=[];
        % TrialListAll(DelI3)=[]; 
     
        ResPoint=find(SLMRes(:,iLaser)==1&sampleN(:,iLaser)<SLMTestParam.TerminalTrialN);  %%Responsive cell needs be further tested        
        DelI4=[];
        for iP=1:length(ResPoint)
            DelItemp1=PointsListAll==ResPoint(iP)&LaserListAll>iLaser;    %%Test current laser
            DelItemp2=PointsListAll==ResPoint(iP)&LaserListAll==iLaser&TrialListAll>SLMTestParam.TerminalTrialN-sampleN(ResPoint(iP),iLaser)+2;    %%but not too many trials within a Tseries                 
            DelI4=union(DelI4,find(DelItemp1+DelItemp2>=1));
        end

        DelTemp1=union(find(DelI1+DelI2+DelI3>=1),DelI4);
        DelI=union(DelI,DelTemp1);
        % DelI4=DelI4>=1;
        % PointsListAll(DelI4)=[];
        % LaserListAll(DelI4)=[];
        % TrialListAll(DelI4)=[];
        
    end
 
    PointsListAll(DelI)=[];
    LaserListAll(DelI)=[];
    TrialListAll(DelI)=[];




    % [LaserListAll,sortI]=sort(LaserListAll);
    % PointsListAll=PointsListAll(sortI);
    % TrialListAll=TrialListAll(sortI);

    [TrialListAll,sortI]=sort(TrialListAll);
    PointsListAll=PointsListAll(sortI);
    LaserListAll=LaserListAll(sortI);


 
    if ZTestNum<=length(PointsListAll)
       Ind=0;
       TestPoints=PointsListAll(Ind+[1:ZTestNum]);
       TestLaserLevels=LaserListAll(Ind+[1:ZTestNum]);
       TestTrial=TrialListAll(Ind+[1:ZTestNum]);

       PointLaserPair=[TestPoints TestLaserLevels];
       % while ~isempty(find(diff(TestPoints)==0))
       %     Ind=Ind+1;
       %     TestPoints=PointsListAll(Ind+[1:ZTestNum]);
       %     TestLaserLevels=LaserListAll(Ind+[1:ZTestNum]);
       % end
    else              %%Need more Points to fullfil the ZTestNum trials in Tseries, using the Points and laser level of responsive cell proved.
       PointLaserPair=unique([PointsListAll LaserListAll],'rows');
       if ~isempty(ResPointLaser)      
          AddPair=ResPointLaser;
       else
          AddP=setdiff([1:N]',PointLaserPair(:,1));
          AddPair=[AddP ones(length(AddP),1)];
       end
       N0=size(PointLaserPair,1);
       N1=size(AddPair,1);


       %%<--------------------------------------------------------------Not
       %%finished yet for adding addtional points while the num of testing
       %%points is not enough
       while size(PointLaserPair,1)<ZTestNum
             if N0>1
                PointLaserPair=[PointLaserPair;AddPair()]
             end

       end
       randperm(N1,min(N1,2));
       %%<--------------------------------------------------------------Not
       %%finished yet for adding addtional points while the num of testing
       %%points is not enough
        if ~isempty(ResPointLaser)      
           
           PointLaserPair=unique([PointsListAll LaserListAll],'rows');
           
           TestPoints=[PointsListAll;ResPointLaser(:,1)];
           TestLaserLevels=[LaserListAll;ResPointLaser(:,2)];
           TestTrial=[TrialListAll;ones(size(ResPointLaser,1),1)];

           TestPoints=repmat(TestPoints,ZTestNum,1);
           TestLaserLevels=repmat(TestLaserLevels,ZTestNum,1);
           TestTrial=repmat(TestTrial,ZTestNum,1);
          
           Ind=0;
           TestPoints=TestPoints(Ind+[1:ZTestNum]);
           TestLaserLevels=TestLaserLevels(Ind+[1:ZTestNum]);
           TestTrial=TestTrial(Ind+[1:ZTestNum]);
     

        end

    end

    figure;hold on;
    plot(TestPoints,'r.-');
    plot(TestLaserLevels,'g^');
    plot(TestTrial,'bo');
    legend({'Point','laser','Trial'})
    % 

end

    


