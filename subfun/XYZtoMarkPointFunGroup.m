function XYZtoMarkPointFunGroup(SavePathAllPoint,Pos3D,IndexNeed,yaml,confSet,Cellstat,FunID)

if nargin<7
   FunID=[];    %%Function ID related to Pos3D(IndexNeed,:).
end

Pos3Dneed=Pos3D(IndexNeed,:);

Repetition=confSet.Repetition;
SpiralSizeUM=confSet.SpiralSizeUM;
SpiralRevolution=confSet.SpiralRevolution;
UncagingLaserPower=confSet.UncagingLaserPower;

SaveName=[SavePathAllPoint 'GPL'];
% SLMIncludedIndFromIscell=IncludedI(tempI);
SLMIncludedIndFromIscell=IndexNeed;
TargetCellstat=Cellstat(SLMIncludedIndFromIscell);
save([SavePathAllPoint 'SLMIncludedIndFromIscell.mat'],'SLMIncludedIndFromIscell','Pos3Dneed','yaml','confSet','Cellstat','TargetCellstat','FunID');



clear Group
GroupCandi=unique(FunID);
for iGroup=1:length(GroupCandi)
    Group(iGroup).Indices=find(FunID==GroupCandi(iGroup));
end
MarkPoints3D_GPLmaker(Pos3Dneed, yaml, true, SpiralSizeUM, SpiralRevolution, SaveName,Group);
% MarkPoints3D_XMLmaker_Points(Pos3Dneed,yaml,true, Repetition, SpiralSizeUM, SpiralRevolution,UncagingLaserPower, SavePathAllPoint);
MarkPoints3D_XMLmaker_Group(Group,Repetition, UncagingLaserPower, SavePathAllPoint);