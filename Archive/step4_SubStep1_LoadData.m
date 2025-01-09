DataFolder=[ProcessFolder 'Data\'];
load([ProcessFolder 'SLMIncludedIndFromIscell.mat'])
[cellIDMap,CellPixCount,MedCenter,cellBoundary,cellIDMapSingleUnit]=Suite2pCellIDMapFromStat(CaData.statCell,[yaml.FOVsize_PX yaml.FOVsize_PY]);
PlaneZ=confSet.ETL+confSet.scan_Z(1);
SLMtargetIDMap=cellIDMapSingleUnit(:,:,SLMIncludedIndFromIscell);
SumDataFolder=[ProcessFolder '\DataSum\'];
mkdir(SumDataFolder)
close all



