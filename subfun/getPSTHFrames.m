function [indexVector, stimulusIDVector, prePostStimuliVector] = getPSTHFrames(InterMPRepetition, PreSLMCal, PostSLMCal,MPFrameJump)
    % GetPSTHFrames - Generates an index vector with PreSLMCal frames before and PostSLMCal frames after each stimulus.
    %
    % Inputs:
    %   InterMPRepetition - Vector containing the cumulative frames where each stimulus occurs.
    %   PreSLMCal - Number of frames before each stimulus.
    %   PostSLMCal - Number of frames after each stimulus.
    %
    % Outputs:
    %   indexVector - A vector containing all the frames that should be captured around each stimulus.
    %   stimulusIDVector - A vector indicating the ID of the stimulus (e.g., 1st, 2nd) for each frame in indexVector.
    %   prePostStimuliVector - A vector indicating whether each frame is before or after the stimulus.

    % Initialize the vectors that will store all indices, stimulus IDs, and pre/post stimuli information
    indexVector = [];
    stimulusIDVector = [];
    prePostStimuliVector = [];

    % Calculate the cumulative frame number for each stimulus
    cumulativeFrames = cumsum(InterMPRepetition);
    % Iterate through each stimulus in cumulativeFrames
    for i = 1:length(cumulativeFrames)-1
        % Get the current stimulus frame
        currentStimulus = cumulativeFrames(i);
        
        if MPFrameJump
       % Calculate frames from PreSLMCal frames before to PostSLMCal frames after the stimulus
        frames = [(currentStimulus - PreSLMCal+1):currentStimulus currentStimulus+2:currentStimulus+PostSLMCal+1];
        else
        % Calculate frames from PreSLMCal frames before to PostSLMCal frames after the stimulus
        frames = (currentStimulus - PreSLMCal):(currentStimulus + PostSLMCal - 1);

        end
        % Add these frames to the index vector
        indexVector = [indexVector, frames];
        
        % Add stimulus ID for each frame
        stimulusIDVector = [stimulusIDVector, repmat(i, 1, length(frames))];
        
        % Add pre/post stimuli information for each frame
        prePostStimuliInfo = [(repmat(-1, 1, PreSLMCal)), (repmat(1, 1, PostSLMCal))];
        prePostStimuliVector = [prePostStimuliVector, prePostStimuliInfo];
    end

    % Check for duplicate frames
    [uniqueIndexVector, uniqueIdx] = unique(indexVector);
    if length(uniqueIndexVector) < length(indexVector)
        warning('Some post-stimuli frames are also pre-stimuli frames for the next stimulus. There are overlapping frames.');
    end

end
