function [FileGenerateInfo,fileList, fileIDs] = getExpInfoFiles_NonMat(folderPath, idRanges)
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
    files = dir(fullfile(folderPath, 'TSeries-*'));
    files = files([files.isdir]==1);
    % Initialize outputs
    allFileList = {files.name}; % Cell array of file names
    allFileIDs = zeros(1, length(allFileList)); % Preallocate for IDs

    for i = 1:length(files)
    % Split the folder name by '-'
       parts = split(allFileList{i}, '-');
    
    % The last part contains the file ID
       allFileIDs(i) = str2num(parts{end});
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

 
    for iFile=1:length(fileList)
        FileGenerateInfo(iFile).FileKey=fileList{iFile};
        FileGenerateInfo(iFile).binFile=[folderPath fileList{iFile} '.bin'];
        FileGenerateInfo(iFile).tifFolder=[folderPath fileList{iFile} '\'];
        FileGenerateInfo(iFile).FileID=fileIDs(iFile);
    end

end
