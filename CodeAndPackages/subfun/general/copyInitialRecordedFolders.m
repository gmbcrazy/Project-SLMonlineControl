function copyInitialRecordedFolders(cellVar, destFolder,varargin)
    % Ensure the destination folder exists
    if nargin==3
       OriginPath=varargin{1};
    else
       OriginPath='';

    end
    
    
    if ~isfolder(destFolder)
        mkdir(destFolder);
    end


    
    % Loop through each folder in the cell variable
    for i = 1:length(cellVar)
        srcFolder = [OriginPath cellVar{i}]; % Get source folder
        
        if isfolder(srcFolder) % Check if source folder exists
            % Extract the folder name correctly
            folderParts = strsplit(srcFolder, filesep);
            folderName = folderParts{end};
            if isempty(folderName)
                folderName = folderParts{end-1};
            end
            % Define destination path
            destPath = fullfile(destFolder, folderName);
            mkdir(destPath);
            % Copy the entire folder recursively
            copyfile(srcFolder, destPath);
            fprintf('Copied folder: %s to %s\n', srcFolder, destPath);
        else
            fprintf('Skipping: %s (Folder does not exist)\n', srcFolder);
        end
    end
end