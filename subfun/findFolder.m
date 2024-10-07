function resultPath = findFolder(rootPath, targetFolderName)
    % findFolder - Recursively searches for a folder in the directory and its subdirectories.
    %
    % Syntax: resultPath = findFolder(rootPath, targetFolderName)
    %
    % Inputs:
    %    rootPath - Path to the root directory where the search should start.
    %    targetFolderName - Name of the folder to search for.
    %
    % Outputs:
    %    resultPath - Path of the found folder. If not found, returns an empty string.

    % Initialize result
    resultPath = '';

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
                % Construct the full path and return it
                resultPath = [fullfile(rootPath, items(i).name) '\'];
                return;
            else
                % Recursively search in subdirectory
                subDirPath = fullfile(rootPath, items(i).name);
                resultPath = findFolder(subDirPath, targetFolderName);
                if ~isempty(resultPath)
                    return;
                end
            end
        end
    end

end
