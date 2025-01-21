function [ExcuteTBL,MatTBL] = FunGroupXMLmatch_BinMat(folderPath, confSet, PSTHparam, Pos3Dneed, idRanges, XMLpattern)
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
%   ExcuteTBL  - Table: Combined output table containing processing results for all files.

% ---------------------------
% STEP 1: Handle optional input for ID ranges
% ---------------------------
if nargin < 6
    idRanges = [];  % If no idRanges is provided, consider all files.
end
MPFrameJump = PSTHparam.MPFrameJump; 
% ---------------------------
% STEP 2: Get the list of files based on ID ranges
% ---------------------------
% Retrieve file names and corresponding IDs that match the given ranges.
[fileList, fileIDs] = getExpInfoFiles(folderPath, idRanges);

% ---------------------------
% STEP 3: Load file generation information
% ---------------------------
% Initialize FileListMAT and load associated metadata.

% Loop through each file and load the `FileGenerateInfo` structure.
for i = 1:length(fileList)
    A = load([folderPath fileList{i}]);  % Load the .mat file
    FileGenerateInfo(i) = A.FileGenerateInfo;  % Extract and store FileGenerateInfo
end

% ---------------------------
% STEP 4: Initialize outputs
% ---------------------------
ExcuteTBL = [];    % Placeholder for combined output tables
MatTBL = [];    % Placeholder for combined output tables
% XMLpattern = 'R(\d+)Laser([\d.]+)GPoint\s?(\d+)';

% ---------------------------
% STEP 5: Process each file
% ---------------------------
for iCount = 1:length(FileGenerateInfo)
    clear OutTBL  % Clear previous output table for new file
    % Extract table for current file using position and configuration settings
    OutTBL = MPSeqFolder_1TargetXNon(FileGenerateInfo(iCount).tifFolder, confSet, Pos3Dneed);
    
    % Assign FileID to the output table for current processing
    OutTBL.FileID = FileGenerateInfo(iCount).FileID * ones(size(OutTBL, 1), 1);
    
    [roundIDs, pointIDs, laserPowers] = XMLPatterExtract(FileGenerateInfo(iCount).xmlFile, XMLpattern);

    TBLmat=[FileGenerateInfo(iCount).FileID * ones(length(roundIDs), 1) roundIDs(:), pointIDs(:), laserPowers(:)];
    MatTBL=[MatTBL;TBLmat];
    % ---------------------------
    % STEP 6: Calculate plane depth and match positions
    % ---------------------------
    % Round plane positions and depths to integers for matching
    planePos = round(Pos3Dneed(OutTBL.Point, 3));  % Extract Z-position (depth)
    planeDepth = round(confSet.ETL + confSet.scan_Z(1));  % Calculate scan plane depths
    
    % Match plane positions to actual depth index
    [~, planeI] = ismember(planePos, planeDepth);
    
    % ---------------------------
    % STEP 7: Generate PSTH (Peri-Stimulus Time Histogram)
    % ---------------------------
    [indexVector, stimulusIDVector, prePostStimuliVector] = getPSTHFrames_MPxmlInfo(...
        OutTBL, PSTHparam.PreSLMCal, PSTHparam.PostSLMCal, PSTHparam.MPFrameJump);
    
    % Calculate PSTH map from the binary file
    [~,InvalidID] = CalMultiPSTHBin(FileGenerateInfo(iCount).binFile, confSet, ...
                              indexVector, stimulusIDVector, prePostStimuliVector);
    OutTBL(InvalidID,:)=[];

    % Append the current table to the aggregated table
    ExcuteTBL = [ExcuteTBL; OutTBL];


    % Clear temporary variables for the current iteration
    clear OutTBL;
end

MatTBL=array2table(MatTBL,'VariableNames',{'FileID','Rounds','Point','UncagingLaserPower',});

end



function [roundIDs, pointIDs, laserPowers] = XMLPatterExtract(MarkPointList, XMLpattern)
    % PatterExtract Extracts round IDs, point IDs, and laser powers from a list of filenames.
    % Inputs:
    %   MarkPointList - A structure array containing file information, typically from a dir() call.
    %   XMLpattern - A regular expression pattern designed to extract specific numerical IDs from the filenames.
    % Outputs:
    %   roundIDs - An array containing numerical round IDs extracted from the file names.
    %   pointIDs - An array containing point identifiers.
    %   laserPowers - An array of laser power values extracted from file names.

    % Initialize arrays to hold extracted data
    roundIDs = zeros(length(MarkPointList), 1);
    laserPowers = zeros(length(MarkPointList), 1);
    pointIDs = zeros(length(MarkPointList), 1);

    % Iterate through each file in the MarkPointList
    for ixml = 1:length(MarkPointList)
        % Extract the file name from the structure
        fileName = MarkPointList{ixml};
        
        % Use regex to parse out the desired data from the file name
        tokens = regexp(fileName, XMLpattern, 'tokens');
        
        % Check if the regex pattern matched and tokens were found
        if ~isempty(tokens)
            % Convert the captured strings to numbers and store them in the respective arrays
            % roundIDs(ixml) = str2double(tokens{1}{1});
            laserPowers(ixml) = str2double(tokens{1}{1});
            pointIDs(ixml) = str2double(tokens{1}{2});
        end
    end

    % Round the extracted numeric values to ensure they are integers
    roundIDs = round(roundIDs);
    pointIDs = round(pointIDs);
end

