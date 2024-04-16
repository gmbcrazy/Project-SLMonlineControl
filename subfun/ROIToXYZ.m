function [SavePath,Pos3D,Pos3DRaw,CaData,CaDataPlane,stat,yaml,confSet]=ROIToXYZ(ConfigFolder)

confSet = ReadYaml([ConfigFolder '\SLMsetting.yml' ]);
LoadPath=confSet.save_path0;
SavePath=[LoadPath 'SingleP\'];
mkdir(SavePath)

ZFocusDepth=confSet.scan_Z;
ZLayer=confSet.ETL+ZFocusDepth;
LoadPath=confSet.save_path0;
Suite2PPath=[LoadPath 'suite2p\'];
xmlFile=dir([LoadPath '*TSeries*-001\TSeries-*-001.xml']);
if isempty(xmlFile)
   disp('Check Path and File Name, No .xml file is detected for recording information')
   SavePath=[];
   Pos3D=[];
   CaData=[];
   stat=[];
   yaml=[];
   confSet=[];
   CaDataPlane=[];
   return;
end
xmlFile=[xmlFile.folder '\' xmlFile.name];
yaml=xml2yaml(xmlFile);
Suite2Temp=[Suite2PPath '\combined\Fall.mat'];
CaData=load(Suite2Temp);
PlaneN=double(CaData.ops.nplanes);
planFolder=dir([Suite2PPath 'plane*']);
clear CaDataPlane;
if length(planFolder)~=PlaneN
   disp('Plane folder number does Not match plane #');

else
    for i=1:PlaneN
    tempFolder=[planFolder(i).folder '\' planFolder(i).name '\Fall.mat'];
    if i==1
       CaDataPlane(i)=load(tempFolder);
    else
       CaDataPlane=concatenateStructs(CaDataPlane, load(tempFolder));
    end

    end
end


Lx=CaDataPlane(1).ops.Lx;
Ly=CaDataPlane(1).ops.Ly;
SLMtargetColor=[0.1 0.8 0.2];

CaDataNew=CaData;
statTemp={};
Shown=[];
CellPlaneID=[];
for i=1:PlaneN
    statTemp=[statTemp CaDataPlane(i).stat];
    CellPlaneID=[CellPlaneID;zeros(length(CaDataPlane(i).stat),1)+i];
    CaData.PlaneMeanImg(:,:,i)=CaDataPlane(i).ops.meanImg;
end
CellPlaneIDRaw=CellPlaneID;

CellPlaneID(CaDataNew.iscell(:,1)==0)=[];
stat=CaDataNew.stat;
stat(CaDataNew.iscell(:,1)==0)=[];
statTemp(CaDataNew.iscell(:,1)==0)=[];


xyPix=[];
for iCell=1:length(statTemp)
    xyPix(iCell,:)=statTemp{iCell}.med(end:-1:1);
end




CaData.CellPlaneID=CellPlaneID;
Pos3D=[];
Zmicro=ZLayer(CellPlaneID);
Pos3D=[xyPix Zmicro(:)];


statRaw=CaDataNew.stat;
xyPixRaw=[];
for iCell=1:length(statRaw)
    xyPixRaw(iCell,:)=statRaw{iCell}.med(end:-1:1);
end

CaData.CellPlaneIDRaw=CellPlaneIDRaw;
Pos3DRaw=[];
Zmicro=ZLayer(CellPlaneIDRaw);
Pos3DRaw=[xyPixRaw Zmicro(:)];



