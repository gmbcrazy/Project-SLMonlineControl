function PSTHmap = CalMultiPSTHBin_Dyn(BinFile, confSet, indexVector, stimulusIDVector, prePostStimuliVector, FrameStep)
    % CalMultiPSTHBin - Calculates the dynamic mean imaging difference (post - pre) for each stimulus.
    %
    % Inputs:
    %   BinFile - Path to the binary file containing imaging data.
    %   confSet - Configuration structure containing imaging parameters.
    %   indexVector - A vector containing all the frames that should be captured around each stimulus.
    %   stimulusIDVector - A vector indicating the ID of the stimulus for each frame in indexVector.
    %   prePostStimuliVector - A vector indicating whether each frame is before or after the stimulus.
    %   FrameStep - Number of frames (window size) to average post-stimulus.
    %
    % Output:
    %   PSTHmap - A map of dimensions [X, Y, nStimuli, WindowNum, nPlane], storing the mean differences.

    % Load imaging data
    ImagingFrame = Suite2pSingleChBin2Frame(BinFile, confSet.SLM_Pixels_Y, confSet.SLM_Pixels_X, length(confSet.ETL), indexVector);

    % Get unique stimulus IDs and count
    uniqueStimuli = unique(stimulusIDVector);
    nStimuli = length(uniqueStimuli);

    % Determine if data is multi-plane
    isMultiPlane = length(confSet.ETL) > 1;
    nPlanes = length(confSet.ETL);

    % Pre-calculate logical masks for pre and post frames
    preMask = (prePostStimuliVector == -1);
    postMask = (prePostStimuliVector == 1);

    % Preallocate PSTHmap
    maxPostFrames = max(sum(stimulusIDVector == uniqueStimuli(1) & postMask));  % Maximum number of post-stimulus frames across stimuli
    WindowNum = ceil(maxPostFrames / FrameStep);  % Maximum number of windows across stimuli

    if isMultiPlane
        PSTHmap = zeros(confSet.SLM_Pixels_Y, confSet.SLM_Pixels_X, nStimuli, WindowNum, nPlanes);
    else
        PSTHmap = zeros(confSet.SLM_Pixels_Y, confSet.SLM_Pixels_X, nStimuli, WindowNum);
    end

    % Loop through each stimulus type
    for i = 1:nStimuli
        stimulusID = uniqueStimuli(i);
        currentMask = (stimulusIDVector == stimulusID);

        % Frame masks for pre- and post-stimulus
        preIndices = currentMask & preMask;
        postIndices = currentMask & postMask;

        % Get the number of post-stimulus frames for this stimulus
        numPostFrames = sum(postIndices);

        if isMultiPlane
            % Multi-plane data
            preImaging = ImagingFrame(:, :, preIndices, :);  % Pre-stimulus frames
            meanPre = mean(preImaging, 3);  % Mean of pre-stimulus frames

            % Loop through dynamic windows
            for w = 1:ceil(numPostFrames / FrameStep)
                frameStart = (w - 1) * FrameStep + 1;
                frameEnd = min(w * FrameStep, numPostFrames);

                postWindowMask = currentMask & postMask;  % Get post-stimulus mask for this stimulus
                postWindowIndices = find(postWindowMask);  % Get the actual frame indices
                postWindowIndices = postWindowIndices(frameStart:frameEnd);  % Select the frame indices for this window

                postImaging = ImagingFrame(:, :, postWindowIndices, :);  % Post-stimulus frames in this window
                meanPost = mean(postImaging, 3);  % Mean of post-stimulus frames in this window

                % Store in PSTHmap
                PSTHmap(:, :, i, w, :) = meanPost - meanPre;
            end
        else
            % Single-plane data
            preImaging = ImagingFrame(:, :, preIndices);  % Pre-stimulus frames
            meanPre = mean(preImaging, 3);  % Mean of pre-stimulus frames

            % Loop through dynamic windows
            for w = 1:ceil(numPostFrames / FrameStep)
                frameStart = (w - 1) * FrameStep + 1;
                frameEnd = min(w * FrameStep, numPostFrames);

                postWindowMask = currentMask & postMask;  % Get post-stimulus mask for this stimulus
                postWindowIndices = find(postWindowMask);  % Get the actual frame indices
                postWindowIndices = postWindowIndices(frameStart:frameEnd);  % Select the frame indices for this window

                postImaging = ImagingFrame(:, :, postWindowIndices);  % Post-stimulus frames in this window
                meanPost = mean(postImaging, 3);  % Mean of post-stimulus frames in this window

                % Store in PSTHmap
                PSTHmap(:, :, i, w) = meanPost - meanPre;
            end
        end
    end
end
