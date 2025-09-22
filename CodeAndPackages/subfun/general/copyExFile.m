% Function: copyExFile
% Purpose: Copy all files from a source folder to a destination folder, skipping files with a certain keyword in their extension
% Author: Lu Zhang
% Date: 8-1-2025

function copyExFile(path, Folder, DestiFolder, ExkeyWord)
    % Handle input arguments
    if nargin < 4
        ExkeyWord = '';
    end
    
    pathFolder = fullfile(path, Folder);
    DestiFolderNew = fullfile(DestiFolder, Folder);
    FileList = dir(pathFolder);
    FileList2 = dir(DestiFolder);
    
    % Check if the destination folder exists, create if not
    CreatFolder = true;
    for i = 1:length(FileList2)
        if FileList2(i).isdir
            if strcmp(Folder, FileList2(i).name)
                CreatFolder = false;
                break;
            end
        end
    end
    if CreatFolder
        mkdir(DestiFolderNew);
    end
    
    % Copy files (excluding ones with ExkeyWord after dot)
    for ii = 1:length(FileList)
        fname = FileList(ii).name;
        if strcmp(fname, '.') || strcmp(fname, '..')
            continue;
        end
        b = strfind(fname, '.');
        if ~isempty(b)
            b = b(end) + 1;
            if b > length(fname)
                continue;
            end
            extStr = fname(b:end);
            if contains(extStr, ExkeyWord)
                continue;
            end
            copyfile(fullfile(pathFolder, fname), DestiFolderNew);
        else
            if ~FileList(ii).isdir
               copyfile(fullfile(pathFolder, fname), DestiFolderNew);
            end
        end
    end
    
    % Recursively copy subfolders
    for ii = 1:length(FileList)
        fname = FileList(ii).name;
        if strcmp(fname, '.') || strcmp(fname, '..')
            continue;
        end
        if FileList(ii).isdir
            % Recursive call
            copyExFile(pathFolder, fname, DestiFolderNew, ExkeyWord);
        end
    end
end
