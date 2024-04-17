function PlotTargetNonTarget(Pos3D,NonTargets,NonTargetsPlane,CaData,NonTargetPath)


NPlane=size(CaData.PlaneMeanImg,3);
GroupXY_Zplane=[Pos3D;NonTargets];
GroupID=[ones(size(Pos3D,1),1);ones(size(NonTargets,1),1)+1];
GroupXY_Zplane(:,3)=[CaData.CellPlaneID(:);NonTargetsPlane(:)];
GroupColor=[0 1 0;1 0 1];
TargetRadius=8;
papersizePX=[0 0 15 15*NPlane];
PlotSLMTargets(CaData.PlaneMeanImg,GroupXY_Zplane,GroupID,GroupColor,TargetRadius);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
saveas(gcf,[NonTargetPath],'fig')
saveas(gcf,[NonTargetPath '.png'],'png')