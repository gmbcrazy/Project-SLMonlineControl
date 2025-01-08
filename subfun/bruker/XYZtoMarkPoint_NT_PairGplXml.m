function XYZtoMarkPoint_NT_PairGplXml(SavePathAllPoint, Pos3D, IndexNeed, NonTargets, IndexNonTargetTrial, yaml, confSet,Cellstat)
% This function processes 3D points for optical stimulation experiments, 
% including selecting specific target and non-target points and generating
% the necessary output files for experimental trials.

% Extract only the needed 3D positions based on IndexNeed.
Pos3DNeed=Pos3D(IndexNeed,:);
% Record the indices of included positions from iscell (indicating selected cells).
SLMIncludedIndFromIscell=IndexNeed;
% Save the relevant indices and position data along with configuration settings to a .mat file.
save([SavePathAllPoint 'SLMIncludedIndFromIscell.mat'],'SLMIncludedIndFromIscell','Pos3DNeed','yaml','confSet','NonTargets','IndexNonTargetTrial','Cellstat');
confSetTemp=confSet;
% Iterate over each trial specified in the configuration settings.
for iTrial=1:confSet.NumTrial
    % Select the non-target points for the current trial using indices provided in IndexNonTargetTrial.
    NonTargetsNeed=NonTargets(IndexNonTargetTrial(:,iTrial),:);
    % Create group points for stimulation, marking both targets and non-targets, and save the result.
%     GroupPoints=MarkPoints3D_GPLmaker_TandNonT(Pos3DNeed, NonTargetsNeed, yaml, true, confSet.SpiralSizeUM, confSet.SpiralRevolution,[SavePathAllPoint 'R' num2str(iTrial)]);
    GroupPoints = MarkPoints3D_GPLmaker_1TAndMultiNonT(Pos3DNeed, NonTargetsNeed, yaml, confSet, [SavePathAllPoint 'R' num2str(iTrial)]);
    % For each stimulation intensity specified, create an XML file for the laser control.
    for iStim =1: length(confSet.UncagingLaserPower)
        confSetTemp.UncagingLaserPower=confSet.UncagingLaserPower(iStim);
        MarkPoints3D_XMLmaker_Group(GroupPoints, confSetTemp,[SavePathAllPoint 'R' num2str(iTrial) 'Laser' num2str(confSet.UncagingLaserPower(iStim))]);
    end
end
