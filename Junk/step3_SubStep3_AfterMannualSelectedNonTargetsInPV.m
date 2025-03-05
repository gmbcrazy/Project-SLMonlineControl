%% Please do mannual correction to exclude disqualified non-targets in PV, after that, exported all selected targets to SelectedFromRaw.gpl file
NonTargetPath=[SavePath  'NonTargetsTest\'];

s=gpl2struct([NonTargetPath 'SelectedFromRaw.gpl']);
% TempTable=struct2table(s.PVGalvoPointList.PVGalvoPoint{1}.Attributes);
for i=1:length(s.PVGalvoPointList.PVGalvoPoint)
    TempStr(i)=s.PVGalvoPointList.PVGalvoPoint{i}.Attributes;
end
NonTargetNeed = convertTableEntries(struct2table(TempStr));
NonTargetNeedInd = NonTargetNeed.Index+1;

NonTargets=NonTargets(NonTargetNeedInd,:);
NonTargetsPlane=NonTargetsPlane(NonTargetNeedInd,:);
SaveNonTargets=[NonTargetPath 'SelectedNonTargets'];
% [NonTargets,NonTargetsPlane]=NonTargetGeneration(SaveNonTargets,NonTargets,NonTargetsPlane,yaml,confSet);
MarkPoints3D_GPLmaker(NonTargets, yaml, confSet,SaveNonTargets,[], 'NonTarget');    %%<< This file is not used in this script, just saved for note, or visulization in PV.

% load([SavePath 'NonTargets\SelectedNonTargets']);

% Visulize all targets and non-targets
figure;
FigSavePath=[NonTargetPath 'TargetNonTarget'];
PlotTargetNonTarget(Pos3D,NonTargets,NonTargetsPlane,CaData,FigSavePath)



% Visulize all targets and non-targets selected for each trial 
for iTrial=1:confSet.NumTrial
    figure;
    IndexNonTargetTrial(:,iTrial)=randperm(size(NonTargets,1),confSet.NumNTperTrial);
    FigSavePath=[NonTargetPath 'Trial' num2str(iTrial)];
    PlotTargetNonTarget(Pos3D,NonTargets(IndexNonTargetTrial(:,iTrial),:),NonTargetsPlane(IndexNonTargetTrial(:,iTrial)),CaData,FigSavePath);
    close all
end


save([NonTargetPath 'NonTarget.mat'],'IndexNonTargetTrial','NonTargets','NonTargetPath','NonTargetsPlane');
save([SavePath 'NonTarget.mat'],'IndexNonTargetTrial','NonTargets','NonTargetPath','NonTargetsPlane');