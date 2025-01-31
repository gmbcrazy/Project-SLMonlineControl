function WorkSpaceFolder = Get_ExpDataFolder(WorkFolder, FolderKeyWord, IncludedList)
    % Function to find an experimental data folder containing specific files
    % 
    % Inputs:
    %   WorkFolder - Root directory to search in
    %   FolderKeyWord - Keyword to filter folders
    %   IncludedList - List of required files/folders within the target folder
    % 
    % Output:
    %   WorkSpaceFolder - Path to the first folder that meets the criteria, empty if none found
    % 
    % Example usage:
    %   WorkSpaceFolder = Get_ExpDataFolder('C:\Data', 'Experiment', {'Data', 'Results.mat'})
    
    WorkSpaceFolder = [];
    resultPaths = findAllFoldersKeyWords(WorkFolder, FolderKeyWord);
    
    for i = 1:length(resultPaths)
        if length(dir(resultPaths{i})) > 2  % Ensures folder is not empty
            item = {dir(resultPaths{i}).name};
            if sum(contains(item, IncludedList)) >= length(IncludedList)
                WorkSpaceFolder = resultPaths{i};
                return;
            end
        end
    end
end
