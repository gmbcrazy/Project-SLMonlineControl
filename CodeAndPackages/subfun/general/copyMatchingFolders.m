function copyMatchingFolders(sourceRoot, sourceFolder, targetRoot, keyword)
% copyMatchingFolders Copies subfolders with a keyword from sourceRoot/sourceFolder to targetRoot/sourceFolder
% Excludes folders containing 'SLMall'
% Ensures destination session folders are created if they don't exist
%
% Usage:
%   copyMatchingFolders('C:\', 'GCamP6S-CamKII', '\\server\path', 'Spon')

    sourcePath = fullfile(sourceRoot, sourceFolder);
    targetPath = fullfile(targetRoot, sourceFolder);
    % Create session folder at destination if it doesn't exist
    if ~exist(targetPath, 'dir')
       mkdir(targetPath);
    end


    % List session folders like SL0886xxxxx
    % sessionFolders = dir(fullfile(sourcePath, 'SL*'));
    % List session folders like SL0886xxxx or L00121xxxxx
    sessionFolders = [dir(fullfile(sourcePath, 'SL*')); dir(fullfile(sourcePath, 'L*'))];


    sessionFolders = sessionFolders([sessionFolders.isdir]);

    for i = 1:length(sessionFolders)
        sessionName = sessionFolders(i).name;
        sessionSource = fullfile(sourcePath, sessionName);
        sessionTarget = fullfile(targetPath, sessionName);

        % Create session folder at destination if it doesn't exist
        if ~exist(sessionTarget, 'dir')
            mkdir(sessionTarget);
        end

        % List subfolders in this session
        subFolders = dir(sessionSource);
        subFolders = subFolders([subFolders.isdir] & ~ismember({subFolders.name}, {'.', '..'}));

        for j = 1:length(subFolders)
            subName = subFolders(j).name;

            if contains(subName, keyword, 'IgnoreCase', true)

                src = fullfile(sessionSource, subName);
                dest = fullfile(sessionTarget, subName);
                fprintf('Copying: %s -> %s\n', src, dest);
                copyfile(src, dest);
            end
        end
    end
end
