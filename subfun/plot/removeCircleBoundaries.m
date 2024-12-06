function cleanedImage = removeCircleBoundaries(image)
    % image: the input medical image
    
    % Convert the image to grayscale if it's not already
    if size(image, 3) == 3
        grayscaleImage = rgb2gray(image);
    else
        grayscaleImage = image;
    end

    % Use Hough Circle Transform to detect circles
    [centers, radii] = imfindcircles(grayscaleImage, [10, 50]);

    % Make a copy of the original image
    cleanedImage = image;

    % Iterate through each detected circle
    for i = 1:size(centers, 1)
        % Get the coordinates and radius of the current circle
        x = centers(i, 1);
        y = centers(i, 2);
        radius = radii(i);

        % Create a circular mask around the circle boundary
        [X, Y] = meshgrid(1:size(image, 2), 1:size(image, 1));
        mask = abs(sqrt((X - x).^2 + (Y - y).^2) - radius) <= 1;

        % Set the intensity values within the circle boundary to a background value (e.g., 0)
        cleanedImage(mask) = 0;
    end
end
