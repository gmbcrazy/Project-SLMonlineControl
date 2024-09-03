function tbl = readGPLFile(filename)
    % readGPLFile - Reads a GPL file and extracts data into a MATLAB table.
    %
    % Syntax: tbl = readGPLFile(filename)
    %
    % Inputs:
    %    filename - The name of the GPL file to read.
    %
    % Outputs:
    %    tbl - A table containing the extracted data.

    % Read the GPL file as text
    fileText = fileread(filename);

    % Parse the XML content from the GPL file
    xDoc = xmlread(fileText);

    % Get the list of PVGalvoPoint elements
    points = xDoc.getElementsByTagName('PVGalvoPoint');

    % Initialize arrays to store data
    X = [];
    Y = [];
    Z = [];
    Name = {};
    Index = [];
    ActivityType = {};
    UncagingLaser = {};
    UncagingLaserPower = [];
    Duration = [];
    IsSpiral = [];
    SpiralSize = [];
    SpiralRevolutions = [];

    % Loop through each point and extract data
    for i = 0:points.getLength-1
        point = points.item(i);
        X = [X; str2double(point.getAttribute('X'))];
        Y = [Y; str2double(point.getAttribute('Y'))];
        Z = [Z; str2double(point.getAttribute('Z'))];
        Name{end+1} = char(point.getAttribute('Name'));
        Index = [Index; str2double(point.getAttribute('Index'))];
        ActivityType{end+1} = char(point.getAttribute('ActivityType'));
        UncagingLaser{end+1} = char(point.getAttribute('UncagingLaser'));
        UncagingLaserPower = [UncagingLaserPower; str2double(point.getAttribute('UncagingLaserPower'))];
        Duration = [Duration; str2double(point.getAttribute('Duration'))];
        IsSpiral = [IsSpiral; strcmp(char(point.getAttribute('IsSpiral')), 'True')];
        SpiralSize = [SpiralSize; str2double(point.getAttribute('SpiralSize'))];
        SpiralRevolutions = [SpiralRevolutions; str2double(point.getAttribute('SpiralRevolutions'))];
    end

    % Create a table from the extracted data
    tbl = table(X, Y, Z, Name', Index, ActivityType', UncagingLaser', ...
                UncagingLaserPower, Duration, IsSpiral, SpiralSize, SpiralRevolutions, ...
                'VariableNames', {'X', 'Y', 'Z', 'Name', 'Index', 'ActivityType', ...
                                  'UncagingLaser', 'UncagingLaserPower', 'Duration', ...
                                  'IsSpiral', 'SpiralSize', 'SpiralRevolutions'});
end
