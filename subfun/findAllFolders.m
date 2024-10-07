function resultPaths = findAllFolders(rootPath, targetFolderName)
    % findAllFolders - Recursively searches for all folders with a specific name 
    % in the directory and its subdirectories.
    %
    % Syntax: resultPaths = findAllFolders(rootPath, targetFolderName)
    %
    % Inputs:
    %    rootPath - Path to the root directory where the search should start.
    %    targetFolderName - Name of the folder to search for.
    %
    % Outputs:
    %    resultPaths - Cell array containing paths of all found folders. 
    %                  If not found, returns an empty cell array.

    % Initialize result
    resultPaths = {};

    % Get list of all items in the current directory
    items = dir(rootPath);
    
    % Loop through each item in the directory
    for i = 1:length(items)
        % Skip '.' and '..' directories
        if strcmp(items(i).name, '.') || strcmp(items(i).name, '..')
            continue;
        end

        % Check if the item is a directory
        if items(i).isdir
            % Check if this is the target folder
            if strcmp(items(i).name, targetFolderName)
                % Construct the full path and add it to the result
                resultPaths{end+1} = [fullfile(rootPath, items(i).name) '\'];
            end
            % Recursively search in subdirectory
            subDirPath = [fullfile(rootPath, items(i).name) '\'];
            subDirResults = findAllFolders(subDirPath, targetFolderName);
            if ~isempty(subDirResults)
                % Append found paths to resultPaths
                resultPaths = [resultPaths, subDirResults];
            end
        end
    end
end
