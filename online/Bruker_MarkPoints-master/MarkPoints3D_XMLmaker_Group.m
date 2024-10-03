function MarkPoints3D_XMLmaker_Group(GroupPoints,Repetition, UncagingLaserPower, SavePath)
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

% get values from settings file


for iGroup=1:length(GroupPoints)

PointIndex = GroupPoints(iGroup).Indices;
if isfield(GroupPoints, 'Name')
    GroupName = GroupPoints(iGroup).Name;
else
    GroupName = ['Group ' num2str(iGroup)];
end
IndicesStr='';
for i=1:length(PointIndex)-1
   IndicesStr=[IndicesStr num2str(PointIndex(i)) ','];
end
IndicesStr=[IndicesStr num2str(PointIndex(end))];


% NumPoints = numel(Xpx);
% for i = 1:NumPoints
  if isfield(GroupPoints, 'PowerWeight')
     CustomLaserPercent=GroupPoints(iGroup).PowerWeight;
     WeightsStr='';
     for i=1:length(PointIndex)-1
         WeightsStr=[WeightsStr num2str(CustomLaserPercent(i)) ','];
     end
     WeightsStr=[WeightsStr num2str(CustomLaserPercent(end))];
     	StimParam = [...
        '  <PVMarkPointElement '...
        'Repetitions="' num2str(Repetition) '" '...
        'UncagingLaser="Satsuma" '...
        'UncagingLaserPower="' num2str(UncagingLaserPower) '" '...
        'CustomLaserPercent="' WeightsStr '" '...
        'TriggerFrequency="' 'None' '" '...
        'TriggerSelection="' 'None' '" '...
        'TriggerCount="' num2str(1) '" '...
        'AsyncSyncFrequency="' 'FirstRepetition' '" '...
        'VoltageOutputCategoryName="None" '...
        'VoltageRecCategoryName="None" '...
        'parameterSet="CurrentSettings" '...
        '>'...

    ];
  else
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

  end



  StimGroup=['    <PVGalvoPointElement '...
        'InitialDelay="10" '...
        'InterPointDelay="5" '...
        'Duration="5" '...
        'SpiralRevolutions="3" '... 
        'Points="' GroupName '" '...
        'Indices="' IndicesStr '" '...
        '/>'...
        ];


% SaveGroupName = ['Laser' num2str(UncagingLaserPower) 'FunGroup' num2str(iGroup)];
% SaveGroupName = ['Laser' num2str(UncagingLaserPower) 'FunGroup' num2str(iGroup)];


% SaveName=[SavePath SaveGroupName];
SaveName=[SavePath GroupName];

SaveName=[SavePath 'Laser' num2str(UncagingLaserPower) GroupName];

if ~strcmpi(SaveName, '')  % if save name provided, save to file
    fid = fopen([SaveName '.xml'], 'w', 'l');
    % fwrite(fid, header, 'char');
    fwrite(fid, '<?xml version="1.0" encoding="utf-8"?>', 'char');
    fprintf(fid,'\n');
    fwrite(fid, '<PVSavedMarkPointSeriesElements Iterations="1" IterationDelay="0.00" CalcFunctMap="False">', 'char');
    fprintf(fid,'\n');
    fwrite(fid, StimParam, 'char');
    fprintf(fid,'\n');
    fwrite(fid, StimGroup, 'char');
    fprintf(fid,'\n');
    fwrite(fid, '  </PVMarkPointElement>', 'char');
    fprintf(fid,'\n');
    fwrite(fid, '</PVSavedMarkPointSeriesElements>', 'char');
    fprintf(fid,'\n');


end
    fclose(fid);

end

