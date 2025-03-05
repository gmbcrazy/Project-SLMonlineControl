%% Load Data
clear all
% ProcessFolder='F:\LuSLMOnlineTest\04222024\SingleP\30PixelFromEdgeExc\';
load('C:\Users\User\Project-SLMonlineControl\subfun\Color\colorMapPN3.mat');
ConfigFolder='C:\Users\User\Project-SLMonlineControl\config\';
ConfigFile='SLMsetting.yml';%<----------------------------------------------------------------------------------Edit, configuration file
[~,~,~,CaData,CaDataPlane,stat,yaml,confSet]=ROIToXYZ(ConfigFolder);
umPerPixel=mean([yaml.umPerlPixelX yaml.umPerlPixelY]);
load('C:\Users\User\Project-SLMonlineControl\subfun\Color\colorMapPN3.mat');
ProcessFolder='F:\LuSLMOnlineTest\SL0242-Ai203\10012024\SingleP\Top13SpeedStimEdgeExc\';%<----------------------Edit, Data folder

step4_SubStep1_LoadData;



%% 
% RandomDelayInterval=[0 1]; %%Random delay is induced after each trial of stimulation.
% PointRepetition=1;  %%Trial Number per each xml MarkPoint stimulation.
PreMarkPointRepetition=25;    %<----------------------------------------------------------------------------------Edit,Frame # before SLM in PV
PostMarkPointRepetition=10;   %<----------------------------------------------------------------------------------Edit,Frame # after SLM in PV
PreSLMCal=15;                 %<----------------------------------------------------------------------------------Edit,Frame # before SLM to calculate baseline map
PostSLMCal=3;                 %<----------------------------------------------------------------------------------Edit,Frame # before SLM to calculate responsive map
% ROIparam.LaserPower=confSet.UncagingLaserPower;  %Laser power to test, using all possible power levels
% ROIparam.LaserPower=confSet.UncagingLaserPower([2 3]);  %<--------------------------------------------------------Edit,It is not necessary to test all possible power levels

step4_SubStep2_PreparingParam;




%% Intiate 1st round test
XMLparam.Laser=1.5;               %<-------------------------------------------------------------------------------Edit, starting laser power to test    
XMLparam.RoundID=1;               %starting round
PointAll=1:size(Pos3Dneed,1);     %All possible MarkPoints for testing
PointsTest=PointAll;              %Initial test Points, this would be updated automatically later
SLMTrialInfo=[];                  %Inital response information, automatically updated after each single trial test
SLMTrialMap=[];                   %Inital response map, automatically updated after each single trial test
clear SLMTable;
SLMTable(:,1)=round(1:size(Pos3Dneed,1));
SLMTable(:,2)=NaN;

FileIDrange=[];                   %<-------------------------------------------------------------------------------Edit, BinFile ID range to calculate SLMresponse
minTrialN=1;
% PointTest=PointAll;
ROIparam.PointsTest=PointsTest;
ROIparam.Clim=[-400 400];
SLMTestParam.TerminalTrialN=5;    %<-------------------------------------------------------------------------------Edit, Trials # to define SLM responsive cells
SLMTestParam.ExcludeTrialN=2;     %<-------------------------------------------------------------------------------Edit, Trials # to define Non-SLM responsive cells
SLMTestParam.AllLaserPower=confSet.UncagingLaserPower;
ROIparam.LaserPower=confSet.UncagingLaserPower;

ROIparam.min_merged_region_size=20;
ROIparam.threshold_percentage=0.25
%% Keep alternating following 2 lines to do all necessary SLM tests and online analysis
tic
step4_SubStep3_InitiateTest
toc

% FileIDrange=[1;400];             %<------------------------------------------------------------------------------Edit, BinFile ID range to calculate SLMresponse
[SLMRes,sampleN]=SLMResponseROIMap(SLMTrialMap,SLMTrialInfo,ROIparam,minTrialN,SumDataFolder,FileIDrange);


ROIparam.LaserPower=confSet.UncagingLaserPower([1 2 3]);  %<-------------------------------------------------------Edit,It is not necessary to test all possible power levels
[UpdateXml, SLMTable, PointsTest, XMLparam, InfoListByLaser]=step3Fun_NextSLMtest(SLMRes,sampleN,ROIparam,XMLparam,SLMTestParam,SLMTable);

PointsTest=[1 2 3 4];
XMLparam.Laser=1.5;
PSTHparam.Plot=1;

%% Save results
save([SumDataFolder 'SLMResponseTable.mat'],'SLMTable','ROIparam','SLMRes','sampleN','SLMTestParam','SLMIncludedIndFromIscell','FunScore');

load([SumDataFolder 'SLMResponseTable.mat'],'SLMTable','ROIparam','SLMRes','sampleN','SLMTestParam','SLMIncludedIndFromIscell','FunScore');


SLMTable([3 5 6 8 15 18 20],2)=1.5
SLMTable([1 2 4 7 11 16 19 25 28 29 30],2)=1.6
SLMTable([17 22 23 24],2)=1.7

TestPower=xmlPower2PVpower(ROIparam.LaserPower);

SMLTablePowerPV=xmlPower2PVpower(SLMTable(:,2));

refPVpower=max(SMLTablePowerPV);
%refPVpower=max(SMLTablePowerPV)+10;  %%make the fixed PV power slightly higher than the maximal power tested, such that we can increase the laser power in real experiment slighly higher than all power we tested.

SMLTablePowerPerc=ceil(PVpower2PVperc(SMLTablePowerPV, refPVpower));

CellPerGroup=10;

FunType=unique(FunScore(:,1));
clear Group
% GroupCandi=unique(FunID);
% for iGroup=1:length(GroupCandi)
%     Group(iGroup).Indices=find(FunID==GroupCandi(iGroup));
%     Group(iGroup).PowerWeight;
% end

FinalPos3D=[];
FinalCellstat={};
FinalFunScore=[];
numFun=length(FunType);


for iGroup=1:numFun
    Group(iGroup).Indices=[1:CellPerGroup]+(iGroup-1)*CellPerGroup;

    FinalFunScore(Group(iGroup).Indices,1)=iGroup;

    Group(iGroup).PowerWeight=ones(CellPerGroup,1);
%     Group(iGroup).PowerWeight(1:length(I1))



    I1=find(FunScore(:,1)==FunType(iGroup)&SLMTable(:,2)>0);
%     Pos3Dtemp=Pos3Dneed(I1,:);
%     cellstattemp=Cellstat(SLMIncludedIndFromIscell(I1));
     if length(I1)>=CellPerGroup  %cells is more than needed
         if iGroup<numFun % functional group cells is mroe than need, choose cells with the highest functional score
            [~,I2]=sort(FunScore(:,iGroup+1),'descend');
            I1=I1(I2(1:CellPerGroup));
         else             % nonfunctional group cells is mroe than need, simple choose the first CellPerGroup cells
            I1=I1(1:CellPerGroup); 
         end
            AddedTargetN=0;
            AddedPos=[];
            AddedCell=cell(1,0);
            AddedFromNonTarget{iGroup}=[];
     else                        %cells is less than needed, added Non-Target as fake cells
          AddedTargetN=CellPerGroup-length(I1);
          AddedPos=NonTargets(1:AddedTargetN,:);
          AddedCell=cell(1,AddedTargetN);
          AddedFromNonTarget{iGroup}=1:AddedTargetN;

     end
      Pos3Dtemp=Pos3Dneed(I1,:);
      cellstattemp=Cellstat(SLMIncludedIndFromIscell(I1));

      Group(iGroup).PowerWeight(1:length(I1))=PVpower2PVweight(xmlPower2PVpower(SLMTable(I1,2)), refPVpower);
%       Group(iGroup).PowerWeight(1:length(I1))=SLMTable(I1,2);

      FinalFunScore(Group(iGroup).Indices(1:length(I1)),2:numFun)=FunScore(I1,2:numFun);
      FinalFunScore(Group(iGroup).Indices(length(I1)+1:CellPerGroup),2:numFun)=nan;
      FinalPos3D=[FinalPos3D;Pos3Dtemp;AddedPos];
      FinalCellstat=[FinalCellstat cellstattemp AddedCell];

end
confSetFinal=confSet;
confSetFinal.UncagingLaserPower=PVpower2xmlPower(refPVpower);
XYZtoMarkPointFunGroup(ProcessFolder,FinalPos3D,Group,yaml,confSetFinal,FinalCellstat)

iFun=1;
find(FunScore(:,1)==FunType(iFun)&SLMTable(:,2)>0)
PointGroup
XYZtoMarkPoint(ProcessFolder,Pos3Dneed,CenterCloseI,yaml,confSet,CaData.statCell);







