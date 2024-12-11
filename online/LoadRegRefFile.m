function [ops, refImg] = LoadRegRefFile(RefFile, FileType, varargin)
% LoadRegRefFile - Load and process reference images for registration.
%
% Inputs:
%   RefFile: Path or folder containing the reference files.
%   FileType: Type of reference file (0: multi-TIFF, 1: Suite2p folder, 2: binary file).
%   varargin: Additional parameters (needed for binary file).
%
% Outputs:
%   ops: Cell array containing registration setup for each image.
%   refImg: 3D matrix of reference images.

% If additional arguments are provided (for binary files)
if nargin == 3
    % Extract binary file parameters
    FileInfo = varargin{1};
    Ly = FileInfo(1);       % Pixel height
    Lx = FileInfo(2);       % Pixel width
    nplane = FileInfo(3);   % Number of planes
    nRepetition = FileInfo(4); % Number of repetitions per plane
end

% Process based on FileType
switch FileType
    case 0 % Multi-TIFF files
        % If no files are provided, open file selection dialog
        if isempty(RefFile)
            [fileName, dirName] = uigetfile('*.tif', 'Select reference image/stack for registration', 'MultiSelect', 'on');
            if ~iscell(fileName)
                fileName = {fileName}; % Ensure fileName is a cell array
            end
        else
            % Get all .tif files in the specified folder
            fileList = dir([RefFile '*.tif']); 
            dirName = fileList.folder; 
            fileName = {fileList.name};
        end

        % Load and process each .tif file
        refImg = [];
        ops = cell(numel(fileName), 1);
        for i = 1:numel(fileName)
            % Read the image and adjust dimensions
            temp = imread([dirName filesep fileName{i}]);
            refImg(:, :, i) = permute(temp, [2 1]); % Adjust for MATLAB indexing
        end

    case 1 % Suite2p folder
        % Get all plane folders in the Suite2p directory
        planeFolder = dir([RefFile 'plane*']);
        nplane = length(planeFolder);
        refImg = [];
        for i = 1:nplane
            % Load mean image for each plane from 'Fall.mat'
            temp = load([planeFolder(i).folder filesep planeFolder(i).name filesep 'Fall.mat'], 'ops');
            refImg(:, :, i) = permute(temp.ops.meanImg,[2 1]); % Adjust for MATLAB indexing
        end

    case 2 % Binary file
        % Convert binary file to frame data and compute mean reference image
        refData = Suite2pSingleChBin2Frame(RefFile, Ly, Lx, nplane, 1:nRepetition);
        refImg = squeeze(mean(refData, 3));

    otherwise
        % Display help message for invalid FileType
        disp('FileType: 0----selected multi tif files or multi tif files folder, 1-------Suite2p folder, 2: bin files');
        return;
end

% Set up registration for each reference image
for i = 1:size(refImg, 3)
    [ops{i, 1}] = setup_registration_phasecorr(refImg(:, :, i));
end
end
