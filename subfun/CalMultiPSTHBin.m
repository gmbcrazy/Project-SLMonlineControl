function [PSTHmap,InvalidStimID] = CalMultiPSTHBin(BinFile, confSet, indexVector, stimulusIDVector, prePostStimuliVector)
    % CalMultiPSTHBin - Calculates the mean imaging difference (post - pre) for each stimulus.
    %
    % Inputs:
    %   BinFile - Path to the binary file containing imaging data.
    %   confSet - Configuration structure containing imaging parameters.
    %   indexVector - A vector containing all the frames that should be captured around each stimulus.
    %   stimulusIDVector - A vector indicating the ID of the stimulus for each frame in indexVector.
    %   prePostStimuliVector - A vector indicating whether each frame is before or after the stimulus.
    %
    % Output:
    %   PSTHmap - A map containing the mean imaging difference (post - pre) for each stimulus.

    % Load imaging data using Suite2pSingleChBin2Frame
    InvalidStimID=[];
    [ImagingFrame,ValidFrame] = Suite2pSingleChBin2Frame(BinFile, confSet.SLM_Pixels_Y, confSet.SLM_Pixels_X, length(confSet.ETL), indexVector);

    InvalidStimID=unique(stimulusIDVector(ValidFrame==0));
    if ~isempty(InvalidStimID)
       DelI=ismember(stimulusIDVector,InvalidStimID);
       indexVector(DelI)=[];
       stimulusIDVector(DelI)=[];
       prePostStimuliVector(DelI)=[];
       if length(size(ImagingFrame))==3      %%1 plane
          ImagingFrame(:,:,DelI)=[];
       elseif length(size(ImagingFrame))==4  %%Multi planes
          ImagingFrame(:,:,DelI,:)=[];
       else
       end
    end

    % Get unique stimulus IDs and count
    uniqueStimuli = unique(stimulusIDVector);
    nStimuli = length(uniqueStimuli);
    
    % Determine if data is multi-plane
    isMultiPlane = length(confSet.ETL) > 1;
    
    % Preallocate PSTHmap
    if isMultiPlane
        PSTHmap = zeros(confSet.SLM_Pixels_Y, confSet.SLM_Pixels_X, nStimuli, length(confSet.ETL));
    else
        PSTHmap = zeros(confSet.SLM_Pixels_Y, confSet.SLM_Pixels_X, nStimuli);
    end

    % Pre-calculate logical masks for pre and post frames
    preMask = (prePostStimuliVector == -1);
    postMask = (prePostStimuliVector == 1);
 
    % Use logical indexing instead of find
    for i = 1:nStimuli
        stimulusID = uniqueStimuli(i);
        currentMask = (stimulusIDVector == stimulusID);

        preIndices = currentMask & preMask;
        postIndices = currentMask & postMask;

        if isMultiPlane
            % Extract and calculate mean for multi-plane data
            preImaging = ImagingFrame(:, :, preIndices, :);
            postImaging = ImagingFrame(:, :, postIndices, :);
            
            meanPre = squeeze(mean(preImaging, 3));
            meanPost = squeeze(mean(postImaging, 3));
            
            PSTHmap(:, :, i, :) = meanPost - meanPre;
        else
            % Extract and calculate mean for single-plane data
            preImaging = ImagingFrame(:, :, preIndices);
            postImaging = ImagingFrame(:, :, postIndices);
            
            meanPre = mean(preImaging, 3);
            meanPost = mean(postImaging, 3);
            
            PSTHmap(:, :, i) = meanPost - meanPre;
        end
    end
end
