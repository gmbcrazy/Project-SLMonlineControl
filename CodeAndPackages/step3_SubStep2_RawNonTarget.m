%% Automatic generate Non-Targets 
NonTargetPath=[SavePath  'NonTargetsTest\'];
mkdir(NonTargetPath)
SaveNonTargets=[NonTargetPath 'Raw'];
BloodVesselTh=0.1;
[NonTargets,NonTargetsPlane]=NonTargetGeneration(SaveNonTargets,Pos3D,CaData.CellPlaneID,yaml,confSet,CaData.PlaneMeanImg,BloodVesselTh);
% confSet.NumNonTarget=20;
% confSet.RadiusAvoidParam=3;
% % [NonTargets,NonTargetsPlane]=NonTargetGeneration(SaveNonTargets,Pos3DRaw,CaData.CellPlaneIDRaw,yaml,confSet);
figure;
FigSavePath=[NonTargetPath 'Raw'];
PlotTargetNonTarget(Pos3D,NonTargets,NonTargetsPlane,CaData,FigSavePath);

