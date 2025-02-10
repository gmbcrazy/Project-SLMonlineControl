function createVideoFromJpeg(folderName, frameRate)
% createVideoFromJpeg Creates a video file from JPEG images in a folder.
%
%   createVideoFromJpeg(folderName, frameRate) reads all .jpeg files in the 
%   specified folder and creates a video file 'Video.avi' using the images 
%   in sorted order with the provided frame rate.
%
%   Inputs:
%       folderName - A string specifying the path to the folder containing the JPEG files.
%       frameRate  - A positive numeric value specifying the video frame rate (frames per second).
%
%   Output:
%       A video file named 'Video.avi' is created in the same folder.
%
%   Note:
%       AVI files can typically be played in Microsoft PowerPoint, provided that 
%       the correct codecs are installed on your system.
%
%   Example:
%       createVideoFromJpeg('C:\MyImages', 30);

    % Check if the folder exists.
    if ~isfolder(folderName)
        error('The folder "%s" does not exist.', folderName);
    end

    % Validate the frameRate input.
    if ~isnumeric(frameRate) || frameRate <= 0
        error('frameRate must be a positive numeric value.');
    end



    % Get a list of all JPEG files in the folder
    jpegFiles = dir(fullfile(folderName, '*.jpg'));
    if isempty(jpegFiles)
        error('No JPEG files found in the folder: %s', folderName);
    end
% Extract filenames
    fileNames = {jpegFiles.name};

    % Natural sorting of filenames
    % Extract numerical part of filenames
    numValues = regexp(fileNames, '\d+', 'match');
    
    % Convert to numerical array (assuming one number per filename)
    numValues = cellfun(@(x) str2double(x{1}), numValues);
    
    % Sort numerically
    [~, sortedIdx] = sort(numValues);
    
    % Reorder jpegFiles based on sorted indices
    jpegFiles = jpegFiles(sortedIdx);
    % Reorder the jpegFiles struct array based on sorted names
    % jpegFiles = jpegFiles(ismember({jpegFiles.name}, fileNames));


    % % Sort the JPEG files alphabetically by name.
    % [~, idx] = sort({jpegFiles.name});
    % jpegFiles = jpegFiles(idx);

    % Define the output video file name and create a VideoWriter object.
    outputFile = fullfile(folderName, 'Video.avi');
    outputVideo = VideoWriter(outputFile);
    outputVideo.FrameRate = frameRate;
    
    % Open the VideoWriter object for writing.
    open(outputVideo);

    % Loop through each JPEG file and write it as a frame in the video.
    for i = 1:length(jpegFiles)
        % Build the full file path.
        filename = fullfile(folderName, jpegFiles(i).name);
        % Read the image.
        img = imread(filename);
        
        % If the image is grayscale, convert it to RGB for a consistent video format.
        if size(img, 3) == 1
            img = repmat(img, [1, 1, 3]);
        end
        
        % Write the current image frame to the video.
        writeVideo(outputVideo, img);
    end

    % Close the video file.
    close(outputVideo);
    
    fprintf('Video created successfully: %s\n', outputFile);
    fprintf('Note: AVI files can generally be played in PowerPoint provided the correct codecs are installed.\n');
end
