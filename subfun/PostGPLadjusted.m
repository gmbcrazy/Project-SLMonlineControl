function PostGPLadjusted(ProcessFolder,IncludeCellFunFilter, FunScore, yaml, confSet,Cellstat,NonTargetPath,NonTargets, IndexNonTargetTrial)

IncludePathOri=[ProcessFolder '\AllIncludedOrigin\'];
IncludePath=[ProcessFolder '\AllIncluded\'];

EditedTable=gpl2Table([IncludePath 'GPL.gpl']);
OriTable=gpl2Table([IncludePathOri 'GPL.gpl']);
% [~,~,I2]=intersect(EditedTable.Name,OriTable.Name);

[~,I2]=ismember(EditedTable.Name,OriTable.Name);



IncludeCellFinal=IncludeCellFunFilter(I2);
if ~isempty(FuncScore)
   FunScore=FunScore(I2,:);
end

SLMIncludedIndFromIscell=IncludeCellFinal;
[XYPosPixel,Z]=gplXYtoPixel(EditedTable,yaml);
Pos3Dneed=[XYPosPixel Z(:,1)];

save([IncludePath 'SLMIncludedIndFromIscell.mat'],'SLMIncludedIndFromIscell','Pos3Dneed','yaml','confSet','NonTargets','IndexNonTargetTrial','Cellstat','FunScore');

% NonTargetPath=

Pos3DFromGPL=[EditedTable.X EditedTable.Y EditedTable.Z];
NonTargetsTable=gpl2Table([NonTargetPath 'SelectedFromRaw.gpl']);
% NonTargetsTable=gpl2Table([NonTargetPath 'Raw.gpl']);

NonTargetsGPL=[NonTargetsTable.X NonTargetsTable.Y NonTargetsTable.Z];

EditedGPLtoMarkPoint_NT_PairGplXml(ProcessFolder, Pos3DFromGPL, SLMIncludedIndFromIscell, NonTargetsGPL, IndexNonTargetTrial, yaml, confSet,Cellstat,FunScore);
