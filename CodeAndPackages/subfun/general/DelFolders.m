function DelFolders(cellVar, varargin) 
    % Ensure the destination folder exists
    if nargin == 2 
       OriginPath = varargin{1};
    else
       OriginPath = '';
    end
    
    % Loop through each folder in the cell variable
    for i = 1:length(cellVar)
        srcFolder = fullfile(OriginPath, cellVar{i}); % Construct full path
        
        if isfolder(srcFolder) % Check if source folder exists
            % Attempt to delete the folder
            try
                rmdir(srcFolder, 's'); % Delete folder and all its contents
                fprintf('Deleted folder: %s\n', srcFolder); 
            catch ME
                fprintf('Failed to delete: %s (Error: %s)\n', srcFolder, ME.message);
            end
        else
            fprintf('Skipping: %s (Folder does not exist)\n', srcFolder);
        end
    end
end
