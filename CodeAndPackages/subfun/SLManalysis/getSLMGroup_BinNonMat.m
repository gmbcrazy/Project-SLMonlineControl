function [PSTHall, OutTBLAll] = getSLMGroup_BinNonMat(folderPath, confSet, PSTHparam, FunGroupPos3D, idRanges)
% GETSLMROI_BINMAT Processes binary .mat files to calculate ROI maps and combines results.
%
% Inputs:
%   folderPath - String: Path to the folder containing ExpInfo files.
%   PSTHparam  - Struct: Parameters for PSTH calculation, including Pre and Post stimulus windows.
%   FunGroupPos3D  -Cell, each Contains 3D position data needed for processing.
%   ROIparam   - Struct: Parameters defining ROI size, like 'NeighbourHfWidthPixel'.
%   idRanges   - Nx2 matrix (optional): Ranges of file IDs to process. If not provided, all files are considered.
%
% Outputs:
%   ROIall     - 3D array: Aggregated ROI heatmaps for all files.
%   OutTBLAll  - Table: Combined output table containing processing results for all files.

% ---------------------------
% STEP 1: Handle optional input for ID ranges
% ---------------------------
if nargin < 5
    idRanges = [];  % If no idRanges is provided, consider all files.
end

% ---------------------------
% STEP 2: Get the list of files based on ID ranges
% ---------------------------
% Retrieve file names and corresponding IDs that match the given ranges.
% [fileList, fileIDs] = getExpInfoFiles(folderPath, idRanges);
[FileGenerateInfo,fileList, fileIDs] = getExpInfoFiles_NonMat(folderPath, idRanges);
% ---------------------------
% STEP 3: Load file generation information
% ---------------------------
% Initialize FileListMAT and load associated metadata.

% Loop through each file and load the `FileGenerateInfo` structure.
% for i = 1:length(fileList)
%     A = load([folderPath fileList{i}]);  % Load the .mat file
%     FileGenerateInfo(i) = A.FileGenerateInfo;  % Extract and store FileGenerateInfo
% end

% ---------------------------
% STEP 4: Initialize outputs
% ---------------------------
ROIall = [];       % Placeholder for aggregated ROI maps
OutTBLAll = [];    % Placeholder for combined output tables
PSTHall =[];
% ---------------------------
% STEP 5: Process each file
% ---------------------------
for iCount = 1:length(FileGenerateInfo)
    clear OutTBL  % Clear previous output table for new file
    FileGenerateInfo(iCount).FileID
    % Extract table for current file using position and configuration settings
    % OutTBL = MPSeqFolder_1TargetXNon(FileGenerateInfo(iCount).tifFolder, confSet, FunGroupPos3D);
    [OutTBL,~]=MPSeqFolder_GroupTargets(FileGenerateInfo(iCount).tifFolder,confSet, FunGroupPos3D);

    % Assign FileID to the output table for current processing
    OutTBL.FileID = FileGenerateInfo(iCount).FileID * ones(size(OutTBL, 1), 1);
    
    % Append the current table to the aggregated table
    OutTBLAll = [OutTBLAll; OutTBL];
    
    % ---------------------------
    % STEP 6: Calculate plane depth and match positions
    % ---------------------------
    % Round plane positions and depths to integers for matching
    % planePos = round(FunGroupPos3D(OutTBL.Point, 3));  % Extract Z-position (depth)
    planeDepth = round(confSet.ETL + confSet.scan_Z(1));  % Calculate scan plane depths
    
    % Match plane positions to actual depth index
    % [~, planeI] = ismember(planePos, planeDepth);
    
    % ---------------------------
    % STEP 7: Generate PSTH (Peri-Stimulus Time Histogram)
    % ---------------------------
    % MPFrameJump = PSTHparam.MPFrameJump;  % Frame jump parameter
    [indexVector, stimulusIDVector, prePostStimuliVector] = getPSTHFrames_MPxmlInfo(...
        OutTBL, PSTHparam.PreSLMCal, PSTHparam.PostSLMCal, PSTHparam.MPFrameJump);
    
    % Calculate PSTH map from the binary file
    PSTHmap = CalMultiPSTHBin(FileGenerateInfo(iCount).binFile, confSet, ...
                              indexVector, stimulusIDVector, prePostStimuliVector);
    
    % Aggregate PSTH maps across files
    PSTHall = cat(3, PSTHall, PSTHmap);
    
    % ---------------------------
    % STEP 8: Calculate ROI Heatmaps
    % ---------------------------
    % for iP = 1:size(OutTBL, 1)  % Loop through all points
    %     Point = OutTBL.Point(iP);  % Extract current point index
    % 
    %     % Generate ROI heatmap for the current point
    %     roiHeat = Center2neighbour(PSTHmap(:, :, iP, planeI(iP)), ...
    %                                FunGroupPos3D(Point, [1, 2]), ...
    %                                ROIparam.NeighbourHfWidthPixel);
    % 
    %     % Smooth the ROI heatmap
    %     roiHeat = SmoothDecDim3(roiHeat, [1 1])';  % Apply smoothing
    % 
    %     % Combine the ROI heatmap into the aggregated ROI array
    %     ROIall = cat(3, ROIall, roiHeat);
    % end
    
    % Clear temporary variables for the current iteration
    clear PSTHmap OutTBL;
end
