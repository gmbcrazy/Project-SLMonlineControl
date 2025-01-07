function MarkPoints3D_XMLmaker_Points(Pos3D,yaml,IsSpiral, Repetition, SpiralSizeUM, SpiralRevolutions, UncagingLaserPower, SavePath)
% Lu Zhang, 02/2024, 
% modifed from MarkPoints_XMLMaker3D.m originally developed by Lloyd Russell 20151119
% Produces XML file of custom galvo positions to be loaded into Prairie View Mark Points

% Noted that this is for multip points stimuli; a group of points could just one single point

% Note:
% unclear why the X uncaging voltages are taken with reference to Galvo X,
% whereas the uncaging Y are taken with reference to Resonant Y...

% Inputs
% ------


% IsSpiral     - string, 'True' for spiral, 'False' for single point
% SpiralSizeUM - size in microns of the sprial
% GroupPoints    - struct array defining groups of points
%                 Example GroupPoints structure
%                 GroupPoints(1).Indices = [1, 2, 3];  % Define the indices
%                 of points belonging to the group, in this case, the 1st, 2nd and 3rd of MarkPoints
%                 GroupPoints(1).Name = 'Group 1';     % Optional: Provide a name for the group
% 
%                 GroupPoints(2).Indices = [4, 5, 6, 7];
%                 GroupPoints(2).Name = 'Group 2';

% SaveName     - provide a save name to save out to file (otherwise returned by function)

% Note pixel coordinates start at 0,0 (0:511 for a 512 image)

% Extract X, Y, and Z coordinates from Pos3D
Xpx=Pos3D(:,1);
Ypx=Pos3D(:,2);
Zmicro=Pos3D(:,3);

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
NumPoints = numel(Xpx);
% Build PointList for individual points

for iStim=1:length(UncagingLaserPower)

	StimParam{iStim} = [...
        '  <PVMarkPointElement '...
        'Repetitions="' num2str(Repetition) '" '...
        'UncagingLaser="Satsuma" '...
        'UncagingLaserPower="' num2str(UncagingLaserPower(iStim)) '" '...
        'TriggerFrequency="' 'None' '" '...
        'TriggerSelection="' 'None' '" '...
        'TriggerCount="' num2str(1) '" '...
        'AsyncSyncFrequency="' 'FirstRepetition' '" '...
        'VoltageOutputCategoryName="None" '...
        'VoltageRecCategoryName="None" '...
        'parameterSet="CurrentSettings" '...
        '>'...

    ];
end

for i = 1:NumPoints
    temp = char(string(IsSpiral));
    IsSpiralStr = [upper(temp(1)) temp(2:end)];
    PointList{i} = ['    <PVGalvoPointElement '...
                    'InitialDelay="10" '...
                    'InterPointDelay="28.3" '...
                    'Duration="5" '...
                    'SpiralRevolutions="' num2str(SpiralRevolutions) '" '...
                    'Points="Point ' num2str(i) '" '...
                    'Indices="' num2str(i) '" '...
                    'X="' num2str(Xv(i)) '" '...
                    'Y="' num2str(Yv(i)) '" '...
                    'IsSpiral="' IsSpiralStr '" '...
                    'SpiralSize="' num2str(SpiralSizeV(end)) '" '...
                    '/>'];
end

for iStim=1:length(UncagingLaserPower)

   for iPoint=1:NumPoints

%       PointIndex = iPoint;
       PointName = ['Laser' num2str(UncagingLaserPower(iStim)) 'Point' num2str(iPoint)];
       SaveName=[SavePath PointName];

     if ~strcmpi(SaveName, '')  % if save name provided, save to file
       fid = fopen([SaveName '.xml'], 'w', 'l');
    % fwrite(fid, header, 'char');
       fwrite(fid, '<?xml version="1.0" encoding="utf-8"?>', 'char');
       fprintf(fid,'\n');
       fwrite(fid, '<PVSavedMarkPointSeriesElements Iterations="1" IterationDelay="0.00" CalcFunctMap="False">', 'char');
       fprintf(fid,'\n');
       fwrite(fid, StimParam{iStim}, 'char');
       fprintf(fid,'\n');
       fwrite(fid,  PointList{iPoint}, 'char');
       fprintf(fid,'\n');
       fwrite(fid, '  </PVMarkPointElement>', 'char');
       fprintf(fid,'\n');
       fwrite(fid, '</PVSavedMarkPointSeriesElements>', 'char');
       fprintf(fid,'\n');
       fclose(fid);

      end

   end
end

