nPlane=length(confSet.ETL)
% DataFolder='F:\LuSLMOnlineTest\04222024\Data\'
DataFolder=[ProcessFolder 'Data\'];
mkdir(DataFolder);
DataLogFolder=[ProcessFolder 'DataLog\'];
SumDataFolder=[ProcessFolder 'DataSum\'];
mkdir(SumDataFolder);

load([ProcessFolder 'SLMIncludedIndFromIscell.mat']);
AllTestPoints3D=Pos3Dneed;
PointAll=1:size(AllTestPoints3D,1);
