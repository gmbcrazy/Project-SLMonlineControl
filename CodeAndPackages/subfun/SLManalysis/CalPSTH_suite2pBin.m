function PSTHall = CalPSTH_suite2pBin(confSet, PSTHparam,TotalMPointInd)
% GETSLMROI_BINMAT Processes binary .mat files to calculate ROI maps and combines results.
%
% Inputs:
%   folderPath - String: Path to the folder containing ExpInfo files.
%   PSTHparam  - Struct: Parameters for PSTH calculation, including Pre and Post stimulus windows.
%   Pos3Dneed  - Numeric array: Contains 3D position data needed for processing.
%   ROIparam   - Struct: Parameters defining ROI size, like 'NeighbourHfWidthPixel'.
%   idRanges   - Nx2 matrix (optional): Ranges of file IDs to process. If not provided, all files are considered.
%
% Outputs:
%   ROIall     - 3D array: Aggregated ROI heatmaps for all files.
%   OutTBLAll  - Table: Combined output table containing processing results for all files.

% ---------------------------
% STEP 1: Handle optional input for ID ranges
% ---------------------------
if nargin < 6
    idRanges = [];  % If no idRanges is provided, consider all files.
end


PSTHall =[];

% ---------------------------
% STEP 5: Process each file
DataFolder=[confSet.save_path0 'suite2p\'];
[indexVector, stimulusIDVector, prePostStimuliVector] = getPSTHFrames_Suite2pBin(TotalMPointInd, PSTHparam.PreSLMCal, PSTHparam.PostSLMCal);

for iplane=1:length(confSet.ETL)
    temppath=[DataFolder 'plane' num2str(iplane-1) '\'];
    tempconfSet=confSet;
    tempconfSet.ETL=0;
    [PSTHmap,InvalidID] = CalMultiPSTHBin([temppath 'data.bin'], tempconfSet, ...
                              indexVector, stimulusIDVector, prePostStimuliVector);
    PSTHall = cat(4, PSTHall, PSTHmap);
    
    % ---------------------------
    % STEP 8: Calculate ROI Heatmaps
    % ---------------------------

    
    % Clear temporary variables for the current iteration
    clear PSTHmap;
end


