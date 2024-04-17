function GroupPoints = MarkPoints3D_GPLmaker_TandNonT(Pos3D, NonTargets, yaml, IsSpiral, SpiralSizeUM, SpiralRevolutions, SaveName)

% Lu Zhang, 02/2024, 
% modifed from MarkPoints_GPLMaker3D.m originally developed by Lloyd Russell 20151119
% This function generates a GPL (Galvo Point List) XML file for custom galvo positions to be loaded into Prairie View Mark Points.
  

% Input:
%   - Pos3D: 3D array of coordinates;
%            Pos3D(:,1) is X coordinates in pixels;
%            Pos3D(:,2) is Y coordinates in pixels;
%            Pos3D(:,3) is Z coordinates in microns;
%   - yaml: structure with recording parameters, typically output of yaml=xml2yaml(xmlFile);
%           xmlFile is the recording-generated XML file from Bruker Prairie View.
%   - IsSpiral: string, 'True' for spiral, 'False' for single point
%   - SpiralSizeUM: size in microns of the spiral
%   - SpiralRevolutions: number of revolutions for spiral points
%   - SaveName: provide a save name to save out to file (otherwise returned by function)
%   - GroupPoints: struct array defining groups of points
%                 Example GroupPoints structure
%                 GroupPoints(1).Indices = [1, 2, 3];  % Define the indices
%                 of points belonging to the group, in this case, the 1st, 2nd and 3rd of MarkPoints
%                 GroupPoints(1).Name = 'Group 1';     % Optional: Provide a name for the group
% 
%                 GroupPoints(2).Indices = [4, 5, 6, 7];
%                 GroupPoints(2).Name = 'Group 2';




% Extract X, Y, and Z coordinates from Pos3D
% Xpx=Pos3D(:,1);
% Ypx=Pos3D(:,2);
% Zmicro=Pos3D(:,3);

Pos3DMix=[Pos3D;NonTargets];
Xpx=Pos3DMix(:,1);
Ypx=Pos3DMix(:,2);
Zmicro=Pos3DMix(:,3);

NumPoint=size(Pos3D,1);
NumNonTarget=size(NonTargets,1);

SingleTargetName={};
for i=1:length(Zmicro)
    if i<=size(Pos3D,1)
       SingleTargetName{i}=['Point ' num2str(i)];
    else
       SingleTargetName{i}=['NonTarget ' num2str(i-size(Pos3D,1))];
    end
end


% Extract relevant parameters from yaml structure
ScanAmp_X = yaml.ScanAmp_X;
ScanAmp_Y = yaml.ScanAmp_Y;
FOVsize_OpticalZoom = yaml.FOVsize_OpticalZoom;
FOVsize_PX = yaml.FOVsize_PX;
FOVsize_UM_1x = yaml.FOVsize_UM;

% Convert full field into imaging FOV
ScanAmp_V_FOV_X = (ScanAmp_X - mean(ScanAmp_X)) + mean(ScanAmp_X);  % Centre, scale, offset
ScanAmp_V_FOV_Y = (ScanAmp_Y - mean(ScanAmp_Y)) + mean(ScanAmp_Y);  % Centre, scale, offset

% Build Look-Up Tables (LUTs)
LUTx = linspace(ScanAmp_V_FOV_X(1), ScanAmp_V_FOV_X(2), FOVsize_PX);
LUTy = linspace(ScanAmp_V_FOV_Y(1), ScanAmp_V_FOV_Y(2), FOVsize_PX);

% Convert pixel coordinates to voltages
Xv = LUTx(Xpx+1);
Yv = LUTy(Ypx+1);
% Xv = LUTx(min(Xpx,yaml.SLM_Pixels_X));
% Yv = LUTy(min(Ypx,yaml.SLM_Pixels_Y));

% Convert spiral size in microns to voltage
SpiralSizeV = SpiralSizeUM * (abs(diff(ScanAmp_V_FOV_X)) / mean(FOVsize_UM_1x));

% Initialize PointList and GroupList
PointList = cell(1, numel(Xpx));

% Build PointList for individual points
for i = 1:numel(Xpx)
    temp = char(string(IsSpiral));
    IsSpiralStr = [upper(temp(1)) temp(2:end)];
    PointList{i} = ['  <PVGalvoPoint '...
                    'X="' num2str(Xv(i)) '" '...
                    'Y="' num2str(Yv(i)) '" '...
                    'Name="' SingleTargetName{i} '" '...
                    'Index="' num2str(i-1) '" '...
                    'ActivityType="MarkPoints" '...
                    'UncagingLaser="Satsuma" '...
                    'UncagingLaserPower="0" '...
                    'Duration="5" '...
                    'IsSpiral="' IsSpiralStr '" '...
                    'SpiralSize="' num2str(SpiralSizeV(end)) '" '...
                    'SpiralRevolutions="' num2str(SpiralRevolutions) '" '...
                    'Z="' num2str(Zmicro(i)) '" '...
                    '/>'];
end

NonTIndex=NumPoint+[1:NumNonTarget];

for iTarget=1:size(Pos3D,1)
    GroupPoints(iTarget).Indices=[iTarget NonTIndex];
    GroupPoints(iTarget).Name=['GPoint ' num2str(iTarget)];
end

GroupList = cell(1, numel(GroupPoints));



% Build GroupList for groups of points
if ~isempty(GroupPoints)
    for i = 1:numel(GroupPoints)
        IndicesStr = '';
        PointIndex = GroupPoints(i).Indices;
        GroupName = GroupPoints(i).Name;
        for j = 1:numel(PointIndex)-1
            IndicesStr = [IndicesStr num2str(PointIndex(j)-1) ','];
        end
        IndicesStr = [IndicesStr num2str(PointIndex(end)-1)];
        GroupList{i} = ['  <PVGalvoPointGroup '...
                        'Indices="' IndicesStr '" '...
                        'Name="' GroupName '" '...
                        'ActivityType="MarkPoints" '...
                        'UncagingLaser="Satsuma" '...
                        'UncagingLaserPower="0" '...
                        'Duration="5" '...
                        'IsSpiral="' IsSpiralStr '" '...
                        'SpiralSize="' num2str(SpiralSizeV(end)) '" '...
                        'SpiralRevolutions="' num2str(SpiralRevolutions) '" '...
                        'Z="' num2str(Zmicro(PointIndex(1))) '" '...
                        '/>'];
    end
end

% Write to file if SaveName provided
if ~strcmpi(SaveName, '')
    fid = fopen([SaveName '.gpl'], 'w', 'l');
    fwrite(fid, '<?xml version="1.0" encoding="utf-8"?>', 'char');
    fprintf(fid,'\n');
    fwrite(fid, '<PVGalvoPointList>', 'char');
    fprintf(fid,'\n');
    for i = 1:numel(Xpx)
        fwrite(fid, PointList{i}, 'char');
        fprintf(fid,'\n');
    end
    if ~isempty(GroupPoints)
        for i = 1:numel(GroupPoints)
            fwrite(fid, GroupList{i}, 'char');
            fprintf(fid,'\n');
        end
    end
    fwrite(fid, '</PVGalvoPointList>', 'char');
    fclose(fid);
end
