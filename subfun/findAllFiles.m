function resultPaths = findAllFiles(rootPath, targetFileName)
    % findAllFiles - Recursively searches for all files with a specific name 
    % in the directory and its subdirectories.
    %
    % Syntax: resultPaths = findAllFiles(rootPath, targetFileName)
    %
    % Inputs:
    %    rootPath - Path to the root directory where the search should start.
    %    targetFileName - Name of the file to search for (including extension).
    %
    % Outputs:
    %    resultPaths - Cell array containing paths of all found files. 
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
            % Recursively search in subdirectory
            subDirPath = fullfile(rootPath, items(i).name);
            subDirResults = findAllFiles(subDirPath, targetFileName);
            if ~isempty(subDirResults)
                % Append found paths to resultPaths
                resultPaths = [resultPaths, subDirResults];
            end
        else
            % Check if the item is the target file
            if strcmp(items(i).name, targetFileName)
                % Construct the full path and add it to the result
                resultPaths{end+1} = fullfile(rootPath, items(i).name);
            end
        end
    end
end
