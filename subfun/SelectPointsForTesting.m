function [TestPoints, TestLaserLevels] = SelectPointsForTesting(SLMRes, sampleN, SLMTestParam, ZTestNum)
    % Get the size of SLMRes and sampleN
    [N, L] = size(SLMRes);
    
    N =10;
    L=3;
    % Generate the PointsListAll and LaserListAll in a vectorized manner
    [iPoints, iTrials, iLasers] = ndgrid(1:N, 1:SLMTestParam.TerminalTrialN, 1:L);
    PointsListAll = iPoints(:);
    LaserListAll = iLasers(:);

    ResPoint=find(sum(SLMRes==1&sampleN>=SLMTestParam.TerminalTrialN,2)==1)
    DelI=ismember(PointsListAll,ResPoint);
    PointsListAll(DelI)=[];
    LaserListAll(DelI)=[];

    for iLaser=1:size(sampleN,2)
        NonResPoint=find(sum(SLMRes(:,iLaser)==0&sampleN(:,iLaser)>=SLMTestParam.ExcludeTrialN,2)==1);
        DelI=ismember(PointsListAll,ResPoint)&LaserListAll=iLaser;
        PointsListAll(DelI)=[];
        LaserListAll(DelI)=[];
    end

    if ZtestNum>=length(PointsListAll)
        TestPoints=PointsListAll(1:ZTestNum);
        TestLaserLevels=LaserListAll(1:ZTestNum);
    end

    

    
end
