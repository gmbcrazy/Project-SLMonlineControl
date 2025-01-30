
LoadPath='\\nih.gov\nimhfileshare\Lab\UFNC\UFNC2\LuZhang\Project3-SLMOnline\LuSLMOnlineTest\01292024\';

caTrialsInd=[1];
% caTrialsInd=[10];

caTrials=XMLread_SLM(LoadPath,caTrialsInd);

xy=[caTrials.AllPoints.X caTrials.AllPoints.Y]


xyPix=floor(xy*512);

Z=[0;0;0;100;100;100];

Pos3D=[xyPix Z];

SpiralSizeUM=15;
SpiralRevolutions=3;
SaveName=['C:\Users\zhangl33\Projects\Project4-SLMonline\LuSLMOnlineTest\01292024\WriteFile\testGPL'];
Group(1).Indices=[1 2];
Group(2).Indices=[1 3 4 6];

MarkPoints_GPLMaker3D_Lu(Pos3D, true, SpiralSizeUM, SpiralRevolutions, SaveName,Group)

SaveName=['C:\Users\zhangl33\Projects\Project4-SLMonline\LuSLMOnlineTest\01292024\WriteFile\testGroup'];

% PointIndex=[1 3 4 6];
Repetiotion=5;
UncagingLaserPower=1.3;

PointIndex=Group(2).Indices;
GroupName='Group 2';
MarkPoints_XMLMaker3DSingleGroup(Pos3D, PointIndex, GroupName,Repetiotion, SpiralSizeUM, UncagingLaserPower,SaveName)


% Lloyd Russell 20151119

thesePoints=Pos3D;
thesePoints(:,4)=1;

        [PhaseMask, TransformedSLMTarget] = SLMPhaseMaskMakerCUDA3D(...
            'Points', thesePoints,...
            'Save', true,...
            'SaveName', [num2str(i, '%03d') '_Z=' num2str(thisZ, '%03d') '_3Dtargets_' timestamp '.tif'],...
            'Do2DTransform', false,...
            'Do3DTransform', false,...
            'AutoAdjustWeights', false);


                [PhaseMask, TransformedSLMTarget] = SLMPhaseMaskMakerCUDA3D(...
            'Points', thesePoints,...
            'Save', false,...
            'Do2DTransform', false,...
            'Do3DTransform', false,...
            'AutoAdjustWeights', false);

[offset_points, zo_position] = ZOBlockAvoider(xyPix)


GPL = MarkPoints_GPLMaker3D([offset_points Z], true, SpiralSizeUM, SpiralRevolutions, SaveName)