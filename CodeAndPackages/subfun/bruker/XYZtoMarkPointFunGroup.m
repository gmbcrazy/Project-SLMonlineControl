function XYZtoMarkPointFunGroup(SavePathAllPoint,Pos3D, Group,yaml,confSet)

if nargin<7
   FunID=[];    %%Function ID related to Pos3D(IndexNeed,:).
end



Repetition=confSet.Repetition;
SpiralSizeUM=confSet.SpiralSizeUM;
SpiralRevolution=confSet.SpiralRevolution;
UncagingLaserPower=confSet.UncagingLaserPower;

SaveName=[SavePathAllPoint 'GPLFunGroup'];
% SLMIncludedIndFromIscell=IncludedI(tempI);
% SLMIncludedIndFromIscell=IndexNeed;
% TargetCellstat=Cellstat(SLMIncludedIndFromIscell);

MarkPoints3D_GPLmaker(Pos3D, yaml, true, SpiralSizeUM, SpiralRevolution, SaveName,Group);
% MarkPoints3D_XMLmaker_Points(Pos3Dneed,yaml,true, Repetition, SpiralSizeUM, SpiralRevolution,UncagingLaserPower, SavePathAllPoint);
MarkPoints3D_XMLmaker_FunGroup(Group,confSet, SavePathAllPoint);

%% Make sham stim xml file
UncagingLaserPower=0.5;
confSetZero=confSet;
confSetZero.UncagingLaserPower=0.5;
MarkPoints3D_XMLmaker_FunGroup(Group,confSetZero, SavePathAllPoint);