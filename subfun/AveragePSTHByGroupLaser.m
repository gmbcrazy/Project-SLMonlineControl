function [averagedPSTHmap,nSample] = AveragePSTHByGroupLaser(XMLTable, PSTHmap,TotalGroupIDs)
    % AveragePSTHByGroupLaser - Averages PSTHmap by different group levels and laser power levels.
    %
    % Inputs:
    %   XMLTable - An n x 2 matrix where the first column contains group IDs and the second column contains laser power levels.
    %   PSTHmap - Imaging data of size [Ly, Lx, Trial, nPlanes].
    %
    % Output:
    %   averagedPSTHmap - A multi-dimensional array of size [Ly, Lx, nPlanes, m], where m is the number of unique group IDs + 1.
    %                     The last element along the 4th dimension represents the average for laserPower == 0.

    % Extract unique group IDs
    groupIDs = unique(XMLTable(:, 1));
    nGroups = length(TotalGroupIDs);
    
    % Determine if data is multi-plane
    isMultiPlane = (length(size(PSTHmap)) == 4);
    [Ly, Lx, ~, nPlanes] = size(PSTHmap);

    % Initialize the averagedPSTHmap array
    averagedPSTHmap = zeros(Ly, Lx, nPlanes, nGroups + 1)+nan;

    % Average for laserPower == 0 (ignoring group IDs)
    matchingTrials = (XMLTable(:, 2) == 0);
    nSampleZero=sum(matchingTrials);
    if isMultiPlane
        selectedPSTHmap = PSTHmap(:, :, matchingTrials, :);
    else
        selectedPSTHmap = PSTHmap(:, :, matchingTrials);
    end
    if ~isempty(selectedPSTHmap)
    meanPSTH = mean(selectedPSTHmap, 3);
    averagedPSTHmap(:, :, :, end) = meanPSTH;
    end

    % Average for laserPower > 0 (by group ID)
    for i = 1:nGroups
        % Get the current group ID
        groupID = TotalGroupIDs(i);
        
        % Find the trials that match the current group ID and laser power > 0
        matchingTrials = (XMLTable(:, 1) == groupID) & (XMLTable(:, 2) > 0);
        nSampleNonZero(i)=sum(matchingTrials);
        % Extract the relevant trials from PSTHmap
        if isMultiPlane
            selectedPSTHmap = PSTHmap(:, :, matchingTrials, :);
        else
            selectedPSTHmap = PSTHmap(:, :, matchingTrials);
        end
        
        % Calculate the mean across the selected trials
        if ~isempty(selectedPSTHmap)
           meanPSTH = mean(selectedPSTHmap, 3);
           averagedPSTHmap(:, :, :, i) = meanPSTH;
        end

    end

    nSample=[nSampleNonZero  nSampleZero];
end
