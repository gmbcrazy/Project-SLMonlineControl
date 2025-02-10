function ImgSmooth = smoothImageSequence(Img, FrameState, windowSize)
% smoothImageSequence Smooths an image sequence using a sliding window mean.
%
%   ImgSmooth = smoothImageSequence(Img, FrameState, windowSize) returns a smoothed 
%   image sequence computed from the input 3D image array Img, a corresponding 1 x n 
%   vector FrameState, and the size of the sliding window (windowSize).
%
%   The smoothing is performed using a sliding window of windowSize consecutive frames,
%   centered on the current frame. Only frames within a contiguous block (i.e. with the 
%   same FrameState value) are averaged together. Frames that do not have enough neighbors 
%   for a full window are left unchanged.
%
%   Inputs:
%     - Img: 3D array of size [height x width x n] (e.g., 512 x 512 x n).
%     - FrameState: 1 x n vector with the state of each frame (e.g., [1 1 1 2 2 2 3 3 3 ...]).
%     - windowSize: Size of the sliding window (must be a positive odd integer).
%
%   Output:
%     - ImgSmooth: 3D array (same size as Img) containing the smoothed images.
%
%   Example:
%     Img = rand(512,512,30);         % Example image stack
%     FrameState = repelem(1:10, 3);    % For instance, 10 states with 3 frames each
%     windowSize = 5;                 % Define the sliding window length
%     ImgSmooth = smoothImageSequence(Img, FrameState, windowSize);

    % Validate that windowSize is a positive odd integer.
    if windowSize < 1 || mod(windowSize, 2) == 0
        error('windowSize must be a positive odd integer.');
    end

    % Get the dimensions of the input image stack.
    [height, width, nFrames] = size(Img);
    
    % Initialize the output with the original images.
    ImgSmooth = Img;
    
    % Determine contiguous segments where FrameState does not change.
    % changeIdx marks the start indices of new segments.
    changeIdx = [1, find(diff(FrameState) ~= 0) + 1, nFrames+1];
    
    % Calculate half window size.
    halfWindow = floor(windowSize / 2);
    
    % Process each contiguous segment independently.
    for seg = 1:length(changeIdx)-1
        segStart = changeIdx(seg);
        segEnd   = changeIdx(seg+1) - 1;
        segIndices = segStart:segEnd;
        L = length(segIndices);
        
        % Only smooth if the segment has at least 'windowSize' frames.
        if L >= windowSize
            % Loop over the frames where a full centered window is available.
            for k = halfWindow+1 : L - halfWindow
                globalIdx = segIndices(k);            % Index in the overall sequence.
                windowIndices = segIndices(k-halfWindow:k+halfWindow);  % Indices for the sliding window.
                % Compute the pixel-wise mean across the specified window.
                ImgSmooth(:,:,globalIdx) = mean(Img(:,:,windowIndices), 3);
            end
        end
        % Frames in segments shorter than windowSize or at the borders remain unchanged.
    end
end
