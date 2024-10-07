function executableList = generateExecutableList(xmlFilesStruct, x, ShamPossibility)
%GENERATEEXECUTABLELIST Generates a list of XML files to be executed.
%   EXECUTABLELIST = GENERATEEXECUTABLELIST(X) returns a cell array containing
%   the names of XML files to be executed X times. The selection is random with
%   Laser0.5 files having ShamPossibility probability and other files having 1-ShamPossibility probability.

    % Get the list of XML files in the current directory
    % xmlFilesStruct = dir('*.xml');
    xmlFiles = {xmlFilesStruct.name};

    % Separate the XML files into categories
    laser05Files = xmlFiles(contains(xmlFiles, 'Laser0.5'));
    laserNonZeroFiles = xmlFiles(~contains(xmlFiles, 'Laser0.5'));
    
    % Probabilities
    probLaser05 = ShamPossibility; % 10% probability for Laser0.5 files
    probLaserNonZero = 1-ShamPossibility; % 90% probability for Laser1.7 files
    
    % Initialize the output list
    executableList = cell(x, 1);
    
    % Generate the executable list
    for i = 1:x
        if rand() < probLaser05
            % Select a random Laser0.5 file
            selectedFile = laser05Files{randi(length(laser05Files))};
        else
            % Select a random Laser1.7 file
            selectedFile = laserNonZeroFiles{randi(length(laserNonZeroFiles))};
        end
        % Add the selected file to the list
        executableList{i} = selectedFile;
    end
end

