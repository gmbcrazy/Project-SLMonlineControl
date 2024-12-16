function [fileList, fileIDs] = getExpInfoFiles(folderPath, idRanges)
% GETEXPINFOFILES Finds 'ExpInfo-*.mat' files in the specified folder and filters them based on ID ranges.
% 
% Inputs:
%   folderPath - String: Path to the folder containing the 'ExpInfo-*.mat' files.
%   idRanges   - Numeric array (optional): Nx2 matrix where each row specifies [minID, maxID] range.
%                If empty or not provided, all files will be considered.
%
% Outputs:
%   fileList - Cell array: Filtered names of 'ExpInfo-*.mat' files.
%   fileIDs  - Numeric array: Extracted file IDs based on the file name.

    % Check if folder exists
    if ~isfolder(folderPath)
        error('The specified folder does not exist.');
    end

    % Get all 'ExpInfo-*.mat' files in the folder
    files = dir(fullfile(folderPath, 'ExpInfo-*.mat'));
    
    % Initialize outputs
    allFileList = {files.name}; % Cell array of file names
    allFileIDs = zeros(1, length(allFileList)); % Preallocate for IDs
    
    % Extract IDs from file names
    for i = 1:length(allFileList)
        % Use regular expression to find numbers after 'ExpInfo-'
        tokens = regexp(allFileList{i}, 'ExpInfo-(\d+)\.mat', 'tokens');
        
        if ~isempty(tokens)
            % Convert the matched number string to a numeric ID
            allFileIDs(i) = str2double(tokens{1}{1});
        else
            allFileIDs(i) = NaN; % In case no valid ID is found
        end
    end

    % Handle optional idRanges input
    if nargin < 2 || isempty(idRanges)
        % If no idRanges are provided, use all files
        validIndices = ~isnan(allFileIDs);
    else
        % Filter files based on idRanges
        validIndices = false(1, length(allFileIDs)); % Initialize logical array
        for r = 1:size(idRanges, 2)
            range = idRanges(:,r);
            validIndices = validIndices | (allFileIDs >= range(1) & allFileIDs <= range(2));
        end
    end

    % Select filtered files and IDs
    fileList = allFileList(validIndices);
    fileIDs = allFileIDs(validIndices);


end
