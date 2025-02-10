function ImgSmooth = ImageSequenceSlidingWin(Img, windowSize)
% smoothImageSequence Smooths an image sequence using a sliding window mean.
%
%   ImgSmooth = smoothImageSequence(Img, windowSize) returns a smoothed image 
%   sequence computed from the input 3D image array Img by averaging over a 
%   sliding window of frames.
%
%   The smoothing is performed using a sliding window of 'windowSize' consecutive 
%   frames centered on the current frame. For frames that do not have enough 
%   neighboring frames for a full window, the original frame value is retained.
%
%   Inputs:
%     - Img: 3D array of size [height x width x nFrames] (e.g., 512 x 512 x nFrames).
%     - windowSize: Size of the sliding window (must be a positive odd integer).
%
%   Output:
%     - ImgSmooth: 3D array (same size as Img) containing the smoothed images.
%
%   Example:
%     Img = rand(512,512,30);  % Example image stack of 30 frames
%     windowSize = 5;          % Define the sliding window length
%     ImgSmooth = smoothImageSequence(Img, windowSize);
%

    % Validate that windowSize is a positive odd integer.
    if windowSize < 1 || mod(windowSize, 2) == 0
        error('windowSize must be a positive odd integer.');
    end

    % Get the dimensions of the input image stack.
    [height, width, nFrames] = size(Img);
    
    % Initialize the output with the original images.
    ImgSmooth = Img;
    
    % Calculate the half-window size.
    halfWindow = floor(windowSize / 2);
    
    % Loop over frames that have enough neighbors for a full window.
    for k = (halfWindow + 1):(nFrames - halfWindow)
        % Define the indices for the current sliding window.
        windowIndices = (k - halfWindow):(k + halfWindow);
        % Compute the mean across the specified window (along the 3rd dimension).
        ImgSmooth(:,:,k) = mean(Img(:,:,windowIndices), 3);
    end
end
