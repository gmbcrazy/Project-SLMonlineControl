function XYZtoMarkPointFunGroup_MultiZ(SavePathAllPoint,Pos3D, Group,yaml,confSet)

if nargin<7
   FunID=[];    %%Function ID related to Pos3D(IndexNeed,:).
end



Repetition=confSet.Repetition;
SpiralSizeUM=confSet.SpiralSizeUM;
SpiralRevolution=confSet.SpiralRevolution;


for iGroup=1:length(Group)

    tempGroup=Group(iGroup);
    tempGroup.Indices=1:length(Group(iGroup).Indices);
    tempGroup.GroupName=['Group 1'];
    % tempGroup.Name=['Group' num2str(iGroup)];
    SaveName=[SavePathAllPoint 'GPLFunGroup' num2str(iGroup)];
% SLMIncludedIndFromIscell=IncludedI(tempI);
% SLMIncludedIndFromIscell=IndexNeed;
% TargetCellstat=Cellstat(SLMIncludedIndFromIscell);

    UncagingLaserPower=confSet.UncagingLaserPower;
    MarkPoints3D_GPLmaker(Pos3D(Group(iGroup).Indices,:), yaml, true, SpiralSizeUM, SpiralRevolution, SaveName,tempGroup);
% MarkPoints3D_XMLmaker_Points(Pos3Dneed,yaml,true, Repetition, SpiralSizeUM, SpiralRevolution,UncagingLaserPower, SavePathAllPoint);
    % MarkPoints3D_XMLmaker_FunGroup(tempGroup,Repetition,UncagingLaserPower, SavePathAllPoint);

    tempGroup.FileName=['Laser' num2str(UncagingLaserPower) 'FunGroup' num2str(iGroup)];
    MarkPoints3D_XMLmaker_Group(tempGroup,Repetition,UncagingLaserPower, SavePathAllPoint);
%% Make sham stim xml file
    UncagingLaserPower=0.5;
    tempGroup.FileName=['Laser' num2str(UncagingLaserPower) 'FunGroup' num2str(iGroup)];
    MarkPoints3D_XMLmaker_Group(tempGroup,Repetition,UncagingLaserPower, SavePathAllPoint);

end