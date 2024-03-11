function GPL = MarkPoints_GPLMaker3D_Lu(Pos3D, IsSpiral, SpiralSizeUM, SpiralRevolutions, SaveName, GroupPoints)
% Lloyd Russell 20151119
% Produces XML file of custom galvo positions to be loaded into Prairie View Mark Points

% Note:
% unclear why the X uncaging voltages are taken with reference to Galvo X,
% whereas the uncaging Y are taken with reference to Resonant Y...

% Inputs
% ------
% Xpx          - 1d array of desired X coordinates (in pixels)
% Ypx          - 1d array of desired Y coordinates (in pixels)
% Pos3D          - 3d array of coordinates; Pos3D(:,1) is X coordinates in
% pixels; Pos3D(:,2) is Y coordinates in pixels; while Pos3D(:,3) is Z
% coordinates in micrometer;

% IsSpiral     - string, 'True' for spiral, 'False' for single point
% SpiralSizeUM - size in microns of the sprial
% SaveName     - provide a save name to save out to file (otherwise returned by function)

% Note pixel coordinates start at 0,0 (0:511 for a 512 image)

% get values from settings file
Xpx=Pos3D(:,1);
Ypx=Pos3D(:,2);
Zmicro=Pos3D(:,3);

yaml = ReadYaml('settings.yml');

% prairie numbers
ScanAmp_X           = yaml.ScanAmp_X;
ScanAmp_Y           = yaml.ScanAmp_Y;
FOVsize_OpticalZoom = yaml.FOVsize_OpticalZoom;
FOVsize_PX          = yaml.FOVsize_PX;
FOVsize_UM_1x       = yaml.FOVsize_UM_1x;

% convert full field into imaging FOV
ScanAmp_V_FOV_X = ((ScanAmp_X - mean(ScanAmp_X)) / FOVsize_OpticalZoom) + mean(ScanAmp_X);  % centre, scale, offset
ScanAmp_V_FOV_Y = ((ScanAmp_Y - mean(ScanAmp_Y)) / FOVsize_OpticalZoom) + mean(ScanAmp_Y);  % centre, scale, offset


% build LUT's
LUTx = linspace(ScanAmp_V_FOV_X(1), ScanAmp_V_FOV_X(2), FOVsize_PX);
LUTy = linspace(ScanAmp_V_FOV_Y(1), ScanAmp_V_FOV_Y(2), FOVsize_PX);

% convert pixel coordinates to voltages
Xv = LUTx(Xpx+1);
Yv = LUTy(Ypx+1);
format long

% convert spiral size in microns to voltage
SpiralSizeV = SpiralSizeUM * (2*ScanAmp_V_FOV_X / (FOVsize_UM_1x/FOVsize_OpticalZoom)); % 0.3V is 10um, 0.6V is 20um

% build the GPL file
% header = [...
%     '<?xml version="1.0" encoding="utf-8"?>\n'...
%     '<PVGalvoPointList>\n'...
%     ];  
NumPoints = numel(Xpx);


% IsSpiral=true;
temp=char(string(IsSpiral));
IsSpiral=[upper(temp(1)) temp(2:end)];
%SpiralSizeV=15;


for i = 1:NumPoints
    % num2str(i)


 PointList{i} = [...
        '  <PVGalvoPoint '...
        'X="' num2str(Xv(i)) '" '...
        'Y="' num2str(Yv(i)) '" '...
        'Name="Point ' num2str(i) '" '...
        'Index="' num2str(i-1) '" '...
        'ActivityType="MarkPoints" '...
        'UncagingLaser="Satsuma" '...
        'UncagingLaserPower="0" '...
        'Duration="10" '...
        'IsSpiral="' IsSpiral '" '...
        'SpiralSize="' num2str(SpiralSizeV(end)) '" '...
        'SpiralRevolutions="' num2str(SpiralRevolutions) '" '...
        'Z="' num2str(Zmicro(i)) '" '...
        '/>'...
    ];

end


if ~isempty(GroupPoints)
    for i = 1:length(GroupPoints)
       IndicesStr='';
      PointIndex=GroupPoints(i).Indices;
      if isfield(GroupPoints,'Name')
         GroupName=GroupPoints(i).Name;
      else
         GroupName=['Group ' num2str(i)];
      end


     for j=1:length(PointIndex)-1
         IndicesStr=[IndicesStr num2str(PointIndex(j)-1) ','];
     end
     IndicesStr=[IndicesStr num2str(PointIndex(end)-1)];

    GroupList{i}=[...
        '  <PVGalvoPointGroup '...
        'Indices="' IndicesStr '" '...
        'Name="' GroupName '" '...
        'ActivityType="MarkPoints" '...
        'UncagingLaser="Satsuma" '...
        'UncagingLaserPower="0" '...
        'Duration="10" '...
        'IsSpiral="' IsSpiral '" '...
        'SpiralSize="' num2str(SpiralSizeV(end)) '" '...
        'SpiralRevolutions="' num2str(SpiralRevolutions) '" '...
        'Z="' num2str(Zmicro(PointIndex(1))) '" '...
        '/>'...
    ];

   end
end



if ~strcmpi(SaveName, '')  % if save name provided, save to file
    fid = fopen([SaveName '.gpl'], 'w', 'l');
    % fwrite(fid, header, 'char');
    fwrite(fid, '<?xml version="1.0" encoding="utf-8"?>', 'char');
    fprintf(fid,'\n');
    fwrite(fid, '<PVGalvoPointList>', 'char');
    fprintf(fid,'\n');

    for i=1:NumPoints
    fwrite(fid, PointList{i}, 'char');
    fprintf(fid,'\n');
    end
    if ~isempty(GroupPoints)
        for i = 1:length(GroupPoints)
           fwrite(fid, GroupList{i}, 'char');
           fprintf(fid,'\n');

        end
    
    end
    fwrite(fid, '</PVGalvoPointList>', 'char');

end
    fclose(fid);

format short
