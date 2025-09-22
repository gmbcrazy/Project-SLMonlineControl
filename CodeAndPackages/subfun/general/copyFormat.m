function copyFormat(path, Folder, DestiFolder, NeedFile, keyWord)
    % copyFormat - Copy specific files (by extension and keyword) recursively.
    %
    % path:         Root directory
    % Folder:       Subfolder name in path to process
    % DestiFolder:  Target folder
    % NeedFile:     File extension to copy (e.g. 'mat', 'txt')
    % keyWord:      Keyword to search in file names (optional)
    
    if nargin < 5
        keyWord = '';
    end

    pathFolder = fullfile(path, Folder);
    DestiFolderNew = fullfile(DestiFolder, Folder);
    FileList = dir(pathFolder);
    FileList2 = dir(DestiFolder);

    % Create destination subfolder if not exist
    CreatFolder = true;
    for i = 1:length(FileList2)
        if FileList2(i).isdir && strcmp(Folder, FileList2(i).name)
            CreatFolder = false;
            break
        end
    end
    if CreatFolder
        mkdir(DestiFolderNew);
    end

    % Copy matching files
    for ii = 1:length(FileList)
        fname = FileList(ii).name;
        if FileList(ii).isdir || strcmp(fname, '.') || strcmp(fname, '..')
            continue
        end
        [~, ~, ext] = fileparts(fname);
        if ~isempty(ext)
            ext = ext(2:end); % remove dot
            if ~isempty(keyWord)
                if strcmp(ext, NeedFile) && contains(fname, keyWord)
                    copyfile(fullfile(pathFolder, fname), DestiFolderNew);
                end
            else
                if strcmp(ext, NeedFile)
                    copyfile(fullfile(pathFolder, fname), DestiFolderNew);
                end
            end
        end
    end

    % Recursively process subfolders
    for ii = 1:length(FileList)
        fname = FileList(ii).name;
        if FileList(ii).isdir && ~strcmp(fname, '.') && ~strcmp(fname, '..')
            copyFormat(pathFolder, fname, DestiFolderNew, NeedFile, keyWord);
        end
    end
end
