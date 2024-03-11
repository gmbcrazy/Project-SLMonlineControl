function GPL = MarkPoints_XMLMaker3DSingleGroup(Pos3D,PointIndex,GroupName,Repetition, SpiralSizeUM, UncagingLaserPower, SaveName)
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


IndicesStr='';
for i=1:length(PointIndex)-1
   IndicesStr=[IndicesStr num2str(PointIndex(i)) ','];
end
IndicesStr=[IndicesStr num2str(PointIndex(end))];

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


% convert spiral size in microns to voltage
SpiralSizeV = SpiralSizeUM * (2*ScanAmp_V_FOV_X / (FOVsize_UM_1x/FOVsize_OpticalZoom)); % 0.3V is 10um, 0.6V is 20um

% build the GPL file
% header = [...
%     '<?xml version="1.0" encoding="utf-8"?>'...
%     '<PVSavedMarkPointSeriesElements>'... 
%     'Iterations="' num2str(p.Results.Iterations) '" '...
%     'IterationDelay="' num2str(p.Results.IterationDelay) '"'...
%     'CalcFunctMap="' "False"...
%     ' >'];  


NumPoints = numel(Xpx);
% for i = 1:NumPoints
	StimParam = [...
        '  <PVMarkPointElement '...
        'Repetitions="' num2str(Repetition) '" '...
        'UncagingLaser="Satsuma" '...
        'UncagingLaserPower="' num2str(UncagingLaserPower) '" '...
        'TriggerFrequency="' 'None' '" '...
        'TriggerSelection="' 'None' '" '...
        'TriggerCount="' num2str(1) '" '...
        'AsyncSyncFrequency="' 'FirstRepetition' '" '...
        'VoltageOutputCategoryName="None" '...
        'VoltageRecCategoryName="None" '...
        'parameterSet="CurrentSettings" '...
        '>'...

    ];

  StimGroup=['    <PVGalvoPointElement '...
        'InitialDelay="10" '...
        'InterPointDelay="5" '...
        'Duration="5" '...
        'SpiralRevolutions="3" '... 
        'Points="' GroupName '" '...
        'Indices="' IndicesStr '" '...
        '/>'...
        ];

% end
% footer = '</PVSavedMarkPointSeriesElements>';
% GPL = [header GroupList footer];

% save the GPL file
% XML = [header GroupList footer];

% save out separate experiment file (may be used in the future?)
% if ~isempty(p.Results.SaveName)
%     fid = fopen([SaveName '.xml'], 'w', 'l');
%     fwrite(fid, XML, 'char');
%     fclose(fid);
% end

if ~strcmpi(SaveName, '')  % if save name provided, save to file
    fid = fopen([SaveName '.xml'], 'w', 'l');
    % fwrite(fid, header, 'char');
    fwrite(fid, '<?xml version="1.0" encoding="utf-8"?>', 'char');
    fprintf(fid,'\n');
    fwrite(fid, '<PVSavedMarkPointSeriesElements Iterations="1" IterationDelay="0.00" CalcFunctMap="False">', 'char');
    fprintf(fid,'\n');
    fwrite(fid, StimParam, 'char');
    fprintf(fid,'\n');
    fwrite(fid,  StimGroup, 'char');
    fprintf(fid,'\n');
    fwrite(fid, '  </PVMarkPointElement>', 'char');
    fprintf(fid,'\n');

    fwrite(fid, '</PVSavedMarkPointSeriesElements>', 'char');
    fprintf(fid,'\n');


end
    fclose(fid);



