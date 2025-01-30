function GPLPoints2XML_SingleMP(GPLfile,PointIndex,SLMconfigYamlFile, SavePath)
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

tbl=gpl2Table(GPLfile);
if isempty(SavePath)
   temp=dir(GPLFile);
   SavePath=[temp.folder '\'];
end


if ~isempty(PointIndex)
   tbl=tbl(PointIndex,:);
end

Yv=tbl.Y;
Xv=tbl.X;
SpiralSizeV=tbl.SpiralSize(1);
SpiralRevolutions=tbl.SpiralRevolutions(1);
IsSpiral=tbl.IsSpiral{1};


if ~isempty(SLMconfigYamlFile)
   if ~isstruct(SLMconfigYamlFile)
   confSet = ReadYaml(SLMconfigYamlFile);
   end
   Repetition=confSet.Repetition;
   UncagingLaserPower=confSet.UncagingLaserPower;
else
   Repetition=8;
   UncagingLaserPower=1.4;
end



% Initialize PointList and GroupList
PointList = cell(1, numel(Xv));
NumPoints = numel(Xv);


MPList=DelBlank(tbl.Name);
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
                    'InterPointDelay="5" '...
                    'Duration="5" '...
                    'SpiralRevolutions="' num2str(SpiralRevolutions) '" '...
                    'Points="' tbl.Name{i} '" '...
                    'Indices="' num2str(tbl.Index(i)+1) '" '...
                    'X="' num2str(Xv(i)) '" '...
                    'Y="' num2str(Yv(i)) '" '...
                    'IsSpiral="' IsSpiralStr '" '...
                    'SpiralSize="' num2str(SpiralSizeV(end)) '" '...
                    '/>'];
end

for iStim=1:length(UncagingLaserPower)

   for iPoint=1:NumPoints

%       PointIndex = iPoint;
       PointName = ['Laser' num2str(UncagingLaserPower(iStim)) MPList{iPoint}];
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


end

function OutputCell=DelBlank(InputCell)

OutputCell=InputCell;
for i=1:length(InputCell)
    I1=findstr(InputCell{i},' ');
    OutputCell{i}(I1)=[];
end
end

