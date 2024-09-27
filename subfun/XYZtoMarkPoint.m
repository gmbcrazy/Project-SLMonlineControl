function XYZtoMarkPoint(SavePathAllPoint,Pos3D,IndexNeed,yaml,confSet,Cellstat,FunScore)

if nargin<7
   FunScore=[];    %%Function score related to Pos3D(IndexNeed,:).
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
save([SavePathAllPoint 'SLMIncludedIndFromIscell.mat'],'SLMIncludedIndFromIscell','Pos3Dneed','yaml','confSet','Cellstat','TargetCellstat','FunScore');



clear Group
Group(1).Indices=[1:length(IndexNeed)];
MarkPoints3D_GPLmaker(Pos3Dneed, yaml, true, SpiralSizeUM, SpiralRevolution, SaveName,Group);
MarkPoints3D_XMLmaker_Points(Pos3Dneed,yaml,true, Repetition, SpiralSizeUM, SpiralRevolution,UncagingLaserPower, SavePathAllPoint);