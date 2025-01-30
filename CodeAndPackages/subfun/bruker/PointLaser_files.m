function OrganizedTable = PointLaser_files(file_names)
    % PointLaser_files - Organizes file names based on laser level, point number, and file ID
    %                    and outputs the result as a table.
    %
    % Syntax: OrganizedTable = PointLaser_files(file_names)
    %
    % Inputs:
    %    file_names - A struct array (from dir) containing file metadata
    %
    % Outputs:
    %    OrganizedTable - A table with columns 'FileNames', 'FileID', 'Laser', and 'Point'
    
    % Regular expression pattern to extract the FileID, laser level, and point number
    pattern = 'TSeries-\d+-\d+-(\d{3})R\dLaser(\d+\.\d+)[a-zA-Z]Point(\d+)\.bin';
    
    % Initialize cell arrays to store data
    fileNamesList = {};
    fileIDList = [];
    laserList = [];
    pointList = [];

    % Loop through each struct in file_names
    for i = 1:length(file_names)
        file_name = file_names(i).name;  % Extract the file name from the struct
        
        % Match the file name against the pattern
        tokens = regexp(file_name, pattern, 'tokens');
        
        if ~isempty(tokens)
            % Extract file ID, laser level, and point number
            file_id = str2double(tokens{1}{1});  % Convert to number
            laser_level = str2double(tokens{1}{2});  % Convert to number
            point_number = str2double(tokens{1}{3});  % Convert to number
            
            % Append data to lists
            fileNamesList{end+1, 1} = file_name;  % File name
            fileIDList(end+1, 1) = file_id;       % File ID as number
            laserList(end+1, 1) = laser_level;    % Laser level as number
            pointList(end+1, 1) = point_number;   % Point number as number
        end
    end

    % Create a table with the organized data
    OrganizedTable = table(fileNamesList, fileIDList, laserList, pointList, ...
                           'VariableNames', {'FileNames', 'FileID', 'Laser', 'Point'});
end
