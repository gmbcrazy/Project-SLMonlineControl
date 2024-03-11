
LoadPath='C:\Users\zhangl33\Projects\Project4-SLMonline\LuSLMOnlineTest\01252024\';


ROIPath='\C:\Users\zhangl33\Projects\Project4-SLMonline\LuSLMOnlineTest\01252024\01252024\1\pointsUpdated.mat';
load(ROIPath);

caTrialsInd=[1];
% caTrialsInd=[10];

caTrials=XMLread_SLM(LoadPath,caTrialsInd);
xmlFile=['C:\Users\zhangl33\Projects\Project4-SLMonline\LuSLMOnlineTest\01252024\TSeries-01252024-001\TSeries-01252024-001.xml'];

yaml=xml2yaml(xmlFile);
xyPix=[points.X(:) points.Y(:)]
% xy
% 
% xyPix=floor(xy*512);
% Z=[0;0;0;100;100;100];

Pos3D=[xyPix points.Zum(:)];

SpiralSizeUM=15;
SpiralRevolutions=3;
SavePath='C:\Users\zhangl33\Projects\Project4-SLMonline\LuSLMOnlineTest\01252024\WriteFile\';
mkdir(SavePath)
SaveName=[SavePath 'GPL'];




Group(1).Indices=[1 2];
Group(2).Indices=[1 3 4 6];

MarkPoints3D_GPLmaker(Pos3D, yaml, true, SpiralSizeUM, SpiralRevolutions, SaveName,Group)

Repetition=5;
UncagingLaserPower=1.275;
MarkPoints3D_XMLmaker_Points(Pos3D,yaml,true, Repetition, SpiralSizeUM, SpiralRevolutions,UncagingLaserPower, SavePath)


MarkPoints3D_XMLmaker_Group(Group,Repetition, SpiralSizeUM, UncagingLaserPower, SavePath)

