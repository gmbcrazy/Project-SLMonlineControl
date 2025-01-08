function GroupPointsAll = MarkPoints3D_GPLMerge_1TAndMultiNonT(Pos3D, NonTargets, yaml, ConfSet, SaveName)

% Lu Zhang, 02/2024, 
% modifed from MarkPoints_GPLMaker3D.m originally developed by Lloyd Russell 20151119
% This function generates a GPL (Galvo Point List) XML file for custom galvo positions to be loaded into Prairie View Mark Points.
  

% Input:
%   - Pos3D: 3D array of coordinates;
%            Pos3D(:,1) is X coordinates in voltage;
%            Pos3D(:,2) is Y coordinates in voltage;
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


SpiralRevolutions=ConfSet.SpiralRevolution;
Duration=ConfSet.Duration;
InterPointDelay=ConfSet.InterPointDelay;
Repetition=ConfSet.Repetition;
IsSpiral=ConfSet.IsSpiral;
SpiralSizeUM=ConfSet.SpiralSizeUM;

% Extract X, Y, and Z coordinates from Pos3D
% Xpx=Pos3D(:,1);
% Ypx=Pos3D(:,2);
% Zmicro=Pos3D(:,3);
NumPoint=size(Pos3D,1);
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
for iPoint=1:NumPoint

Pos3DMix=[Pos3D(iPoint,:);NonTargets];

Zmicro=Pos3DMix(:,3);

NumNonTarget=size(NonTargets,1);

SingleTargetName={};

% SingleTargetName{1}=['Point ' num2str(iPoint)];
% for i=2:length(Zmicro)
%     SingleTargetName{i}=['NonTarget ' num2str(i-1)];
% end
% 
% For multi Zseries, it is good to have standarized naming
for i=1:length(Zmicro)
    SingleTargetName{i}=['Point ' num2str(i)];
end


Xv=Pos3DMix(:,1);
Yv=Pos3DMix(:,2);


% Convert spiral size in microns to voltage
SpiralSizeV = SpiralSizeUM * (abs(diff(ScanAmp_V_FOV_X)) / mean(FOVsize_UM_1x));

% Initialize PointList and GroupList
PointList = cell(1, numel(Xv));

% Build PointList for individual points
for i = 1:numel(Xv)
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

NonTIndex=1+[1:NumNonTarget];


clear GroupPoints


%     GroupPoints(1).Indices=[1 NonTIndex];
%     GroupPoints(1).Name=['Point' num2str(iPoint)];
%     GroupPoints(1).FileName=['GPoint' num2str(iPoint)];

% For multi Zseries, it is good to have standarized naming
    clear GroupPoints
    GroupPoints(1).Indices=[1 NonTIndex];
    GroupPoints(1).Name=['Group ' num2str(1)];
    GroupPoints(1).FileName=['GPoint' num2str(iPoint)];


GroupList = cell(1, numel(GroupPoints));



% Build GroupList for groups of points
if ~isempty(GroupPoints)
    for i = 1:numel(GroupPoints)
        IndicesStr = '';
        PointIndex = GroupPoints(i).Indices;
        GroupName = GroupPoints(i).Name;
%         GroupName = GroupPoints(i).FileName;
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
    fid = fopen([SaveName GroupPoints(1).FileName '.gpl'], 'w', 'l');
    fwrite(fid, '<?xml version="1.0" encoding="utf-8"?>', 'char');
    fprintf(fid,'\n');
    fwrite(fid, '<PVGalvoPointList>', 'char');
    fprintf(fid,'\n');
    for i = 1:numel(Xv)
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

GroupPointsAll(iPoint)=GroupPoints;

end
