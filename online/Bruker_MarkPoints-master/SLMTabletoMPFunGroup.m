function SLMTabletoMPFunGroup(SavePathAllPoint,Pos3D,SLMTable,yaml,confSet,Cellstat,FunID)

if nargin<7
   FunID=[];    %%Function ID related to Pos3D(IndexNeed,:).
end
IndexNeed=find(SLMTable(:,2)>1);   %%SLM responseive cells is included only


Pos3Dneed=Pos3D(IndexNeed,:);

Repetition=confSet.Repetition;
SpiralSizeUM=confSet.SpiralSizeUM;
SpiralRevolution=confSet.SpiralRevolution;
UncagingLaserPower=confSet.UncagingLaserPower;

SaveName=[SavePathAllPoint 'GPL'];

SLMIncludedIndFromSLMTable=IndexNeed;
TargetCellstat=Cellstat(SLMIncludedIndFromSLMTable);
save([SavePathAllPoint 'SLMIncludedIndFromSLMTable.mat'],'SLMTable','SLMIncludedIndFromSLMTable','Pos3Dneed','yaml','confSet','TargetCellstat','FunID');



clear Group
GroupCandi=unique(FunID);
for iGroup=1:length(GroupCandi)
    Group(iGroup).Indices=find(FunID==GroupCandi(iGroup));
    GroupPoint(iGroup).PowerWeight;
end
MarkPoints3D_GPLmaker(Pos3Dneed, yaml, true, SpiralSizeUM, SpiralRevolution, SaveName,Group);
% MarkPoints3D_XMLmaker_Points(Pos3Dneed,yaml,true, Repetition, SpiralSizeUM, SpiralRevolution,UncagingLaserPower, SavePathAllPoint);
MarkPoints3D_XMLmaker_Group(Group,Repetition, UncagingLaserPower, SavePathAllPoint);