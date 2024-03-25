%% 
LoadPath='F:\LuSLMOnlineTest\03222024\';


ZFocusDepth=147.17;
ZLayer=[-50 0 50]+ZFocusDepth;

Suite2PPath=[LoadPath 'suite2p\']
mkdir(Suite2PPath)
xmlFile=dir([LoadPath '*TSeries*-001\TSeries-*-001.xml'])

xmlFile=[xmlFile.folder '\' xmlFile.name];
yaml=xml2yaml(xmlFile);

% xmlFile=['F:\LuSLMOnlineTest\02272024\TSeries-02272024-1114-001\TSeries-02272024-1114-001.xml'];

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
end
CellPlaneID(CaDataNew.iscell(:,1)==0)=[];
stat=CaDataNew.stat;
stat(CaDataNew.iscell(:,1)==0)=[];
statTemp(CaDataNew.iscell(:,1)==0)=[];


xyPix=[];
for iCell=1:length(statTemp)
    xyPix(iCell,:)=statTemp{iCell}.med(end:-1:1);
end











% ROIPath='\C:\Users\zhangl33\Projects\Project4-SLMonline\LuSLMOnlineTest\01252024\01252024\1\pointsUpdated.mat';
% load(ROIPath);

% caTrialsInd=[1];
% caTrialsInd=[10];

% caTrials=XMLread_SLM(LoadPath,caTrialsInd);



Pos3D=[];
% PointPath='F:\LuSLMOnlineTest\02082024\02082024\1\*POINT*.mat';
% A=dir(PointPath)




Zmicro=ZLayer(CellPlaneID);
% 

Pos3D=[xyPix Zmicro(:)] 

SLMrangePix=50;
XYrange=[SLMrangePix;512-SLMrangePix]  %%Cell locates close to edge of the view, were not considered as SML targets.
OutRange=find(xyPix(:,1)<XYrange(1)|xyPix(:,2)<XYrange(1)|xyPix(:,1)>XYrange(2)|xyPix(:,2)>XYrange(2))
IncludedI=setdiff(1:length(Zmicro),OutRange);
Pos3D(OutRange,:)=[];


tempI=1:1:size(Pos3D,1);


%  2 4 7 13 14 15 18 23 27

tempI=[4 8 13 14 18 26 28]

Pos3Dneed=Pos3D(tempI,:)
SpiralSizeUM=15;
SpiralRevolutions=3;
SavePath=[LoadPath '\SingleP\'];
mkdir(SavePath)
SaveName=[SavePath 'GPL'];
SLMIncludedIndFromIscell=IncludedI(tempI);
save([SavePath 'SLMIncludedIndFromIscell.mat'],'SLMIncludedIndFromIscell')







%%

clear Group
Group(1).Indices=[1:length(tempI)];
% Group(2).Indices=[10:19];
% Group(3).Indices=[1:5]+45;
% Group(4).Indices=[1:3:60];

MarkPoints3D_GPLmaker(Pos3Dneed, yaml, true, SpiralSizeUM, SpiralRevolutions, SaveName,Group)

Repetition=5;
UncagingLaserPower=[1.12 1.14];
MarkPoints3D_XMLmaker_Points(Pos3Dneed,yaml,true, Repetition, SpiralSizeUM, SpiralRevolutions,UncagingLaserPower, SavePath)


% MarkPoints3D_XMLmaker_Group(Group,Repetition, SpiralSizeUM, UncagingLaserPower, SavePath)

