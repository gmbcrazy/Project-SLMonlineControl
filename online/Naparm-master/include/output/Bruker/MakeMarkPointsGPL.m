function MakeMarkPointsGPL(handles)
GalvoPositions    = handles.GalvoPositions;
X                 = GalvoPositions(:,1);
Y                 = GalvoPositions(:,2);
Group_info=handles.points;%%added by Hari
SpiralRevolutions = str2double(handles.SpiralRevolutions_Edit.String);
SpiralDiameterUm  = str2double(handles.SpiralDiameter_Edit.String);
IsSpiral          = 'True';
SaveName          = [handles.data.NaparmDirectory filesep handles.data.ExperimentIdentifier];

MarkPoints_GPLMaker(X, Y, IsSpiral, SpiralDiameterUm, SpiralRevolutions, SaveName,Group_info);

