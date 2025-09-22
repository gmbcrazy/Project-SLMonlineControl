function [Pos3D,Pos3DRaw,CaData,CaDataPlane,stat]=Extract_Suite2p(confSet)


LoadPath=confSet.save_path0;
% SavePath=[LoadPath 'SingleP\'];
% mkdir(SavePath)

ZFocusDepth=confSet.scan_Z;
ZLayer=confSet.ETL+ZFocusDepth;
LoadPath=confSet.save_path0;
Suite2PPath=[LoadPath 'suite2p\'];

Suite2Temp=[Suite2PPath '\combined\Fall.mat'];
CaData=load(Suite2Temp);
PlaneN=double(CaData.ops.nplanes);
planeFolder=dir([Suite2PPath 'plane*']);
clear CaDataPlane;
if length(planeFolder)~=PlaneN
   disp('Plane folder number does Not match plane #');

else
    for i=1:PlaneN
    tempFolder=[planeFolder(i).folder '\' planeFolder(i).name '\Fall.mat'];
    if i==1
       CaDataPlane(i)=load(tempFolder);
    else
       CaDataPlane=concatenateStructs(CaDataPlane, load(tempFolder));
    end

    end
end
CaData.F=double(CaData.F);
CaData.spks=double(CaData.spks);
CaData.Fneu=double(CaData.Fneu);

CaData.dF=CaData.F - CaData.ops.neucoeff * CaData.Fneu;


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
    CaData.PlaneMeanImg(:,:,i)=double(CaDataPlane(i).ops.meanImg);
    CaData.PlaneMeanImgE(:,:,i)=double(CaDataPlane(i).ops.meanImgE);

end

statRaw=statTemp;

CellPlaneIDRaw=CellPlaneID;

CellPlaneID(CaDataNew.iscell(:,1)==0)=[];
stat=CaDataNew.stat;
stat(CaDataNew.iscell(:,1)==0)=[];
statTemp(CaDataNew.iscell(:,1)==0)=[];


xyPix=[];
for iCell=1:length(statTemp)
    xyPix(iCell,:)=statTemp{iCell}.med(end:-1:1);
end



CaData.statCell=statTemp;

CaData.CellPlaneID=CellPlaneID;
Pos3D=[];
Zmicro=ZLayer(CellPlaneID);
Pos3D=[xyPix Zmicro(:)];


xyPixRaw=[];
for iCell=1:length(statRaw)
    xyPixRaw(iCell,:)=statRaw{iCell}.med(end:-1:1);
end

CaData.CellPlaneIDRaw=CellPlaneIDRaw;
Pos3DRaw=[];
Zmicro=ZLayer(CellPlaneIDRaw);
Pos3DRaw=[xyPixRaw Zmicro(:)];



