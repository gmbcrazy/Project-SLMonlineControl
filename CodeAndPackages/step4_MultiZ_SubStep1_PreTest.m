[PointLaserPair,ResPointLaser] = SelectPointsForTesting_v4(SLMRes, sampleN, SLMTestParam, PowerTestPVPar.Ziteration-1);
disp('Next Point and Power Level Index (ROI.LaserPower) for Test:')
disp(PointLaserPair')

SLMTable(ResPointLaser(:,1),2)=ROIparam.LaserPower(ResPointLaser(:,2))';

if isempty(PointLaserPair)
   disp('Terminate the power test')
else
   XMLparam.PointList=PointLaserPair(:,1);                   %%nP, NumOfTestedPoints, nP + 1 = NumOfZseries 
%    ROIparam.PointsTest=XMLparam.PointList;
   XMLparam.Laser=ROIparam.LaserPower(PointLaserPair(:,2));  %Noted that, laser test levels is dependent on ROIparam.LaserPower, not SLMTestParam.AllLaserPower
   XMLparam.RoundID=randperm(XMLparam.TotalRounds,1);    %next round of Non-targets were chosen.
   PSTHparam.PointsTest=ROIparam.PointsTest;
   PSTHparam.LaserPower=ROIparam.LaserPower;

   PVpower=xmlPower2PVpower(XMLparam.Laser(1));
   disp(['Change Satusma Power to ' num2str(PVpower)]);

end

if sum(abs([length(XMLparam.Laser) length(XMLparam.PointList)]+1-PowerTestPVPar.Ziteration))~=0
   disp('Check whether # Point, Laser levels and Zseries match')
end