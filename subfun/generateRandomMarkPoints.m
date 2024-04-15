function markPoints = generateRandomMarkPoints(numPoints, numPlanes, planeSize)
    % Initialize the array to store MarkPoints
    markPoints = zeros(numPoints, 3);
    
    % Randomly generate MarkPoints
    for i = 1:numPoints
        x = randi(planeSize(2));  % Random x-coordinate within the plane width
        y = randi(planeSize(1));  % Random y-coordinate within the plane height
        z = randi(numPlanes);     % Random z-coordinate (plane number)
        
        % Store the generated point
        markPoints(i, :) = [x, y, z];
    end
end


