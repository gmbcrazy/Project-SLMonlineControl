function deleteFormat(pathFolder, DeletFile, ExcludeKeyWord)
% Recursively delete files of a certain format in all subfolders,
% but skip files whose names include ExcludeKeyWord.
%
% Usage:
%   deleteFormat(pathFolder, DeletFile, ExcludeKeyWord)
% Example:
%   deleteFormat('C:/Data', 'mat', 'keep')
%   % Deletes all .mat files except those with 'keep' in the name

    if nargin < 3
        ExcludeKeyWord = '';
    end

    % List all files/folders
    FileList = dir(pathFolder);

    % Delete matching files in current folder (excluding ExcludeKeyWord)
    for i = 1:length(FileList)
        name = FileList(i).name;
        filePath = fullfile(pathFolder, name);

        if ~FileList(i).isdir
            [~, ~, ext] = fileparts(name);
            if length(ext) > 1 && strcmpi(ext(2:end), DeletFile)
                if isempty(ExcludeKeyWord) || ~contains(name, ExcludeKeyWord)
                    delete(filePath);
                end
            end
        end
    end

    % Recurse into subfolders
    for i = 1:length(FileList)
        name = FileList(i).name;
        if FileList(i).isdir && ~strcmp(name, '.') && ~strcmp(name, '..')
            deleteFormat(fullfile(pathFolder, name), DeletFile, ExcludeKeyWord);
        end
    end
end
