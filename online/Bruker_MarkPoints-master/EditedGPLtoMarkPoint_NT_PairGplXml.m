function EditedGPLtoMarkPoint_NT_PairGplXml(SavePathAllPoint, Pos3DFromGPL, SLMIncludedIndFromIscell, NonTargetsGPL, IndexNonTargetTrial, yaml, confSet,Cellstat,FunScore)
% This function processes 3D points for optical stimulation experiments, 
% including selecting specific target and non-target points and generating
% the necessary output files for experimental trials.
if nargin<9
   FunScore=[];    %%Function score related to Pos3D(IndexNeed,:).
end
% Extract only the needed 3D positions based on IndexNeed.
[XYPosPixel,Z]=gplXYPostoPixel(Pos3DFromGPL,yaml);
Pos3DNeed=[XYPosPixel Z(:,1)];
[NonTargetsPosPixel,ZNonTargets]=gplXYPostoPixel(NonTargetsGPL,yaml);
NonTargets=[NonTargetsPosPixel ZNonTargets(:,1)];

% Pos3DNeed=Pos3D(IndexNeed,:);
% Record the indices of included positions from iscell (indicating selected cells).
% Save the relevant indices and position data along with configuration settings to a .mat file.
save([SavePathAllPoint 'SLMIncludedIndFromIscell.mat'],'SLMIncludedIndFromIscell','Pos3DNeed','yaml','confSet','NonTargets','IndexNonTargetTrial','Cellstat','FunScore');

% Iterate over each trial specified in the configuration settings.
for iTrial=1:confSet.NumTrial
    % Select the non-target points for the current trial using indices provided in IndexNonTargetTrial.
    NonTargetsNeed=NonTargetsGPL(IndexNonTargetTrial(:,iTrial),:);
    % Create group points for stimulation, marking both targets and non-targets, and save the result.
%     GroupPoints=MarkPoints3D_GPLmaker_TandNonT(Pos3DNeed, NonTargetsNeed, yaml, true, confSet.SpiralSizeUM, confSet.SpiralRevolution,[SavePathAllPoint 'R' num2str(iTrial)]);
%     GroupPoints = MarkPoints3D_GPLmaker_1TAndMultiNonT(Pos3DNeed, NonTargetsNeed, yaml, true, confSet.SpiralSizeUM, confSet.SpiralRevolution, [SavePathAllPoint 'R' num2str(iTrial)]);
    GroupPoints = MarkPoints3D_GPLMerge_1TAndMultiNonT(Pos3DFromGPL, NonTargetsNeed, yaml, true, confSet.SpiralSizeUM, confSet.SpiralRevolution, [SavePathAllPoint 'R' num2str(iTrial)]);
    % For each stimulation intensity specified, create an XML file for the laser control.
    for iStim =1: length(confSet.UncagingLaserPower)
        MarkPoints3D_XMLmaker_Group(GroupPoints, confSet.Repetition, confSet.UncagingLaserPower(iStim),[SavePathAllPoint 'R' num2str(iTrial) 'Laser' num2str(confSet.UncagingLaserPower(iStim))]);
    end
end
