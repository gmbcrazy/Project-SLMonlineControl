function [Neighbourhood,  nonNeighbourhood] = MarkPoint2Neighbourhood(MarkPoints, radius, numPlanes,planeSize)
    % Define the size of each plane
    % planeSize = [Ly, Lx];
    
    % Initialize the 3D matrix with 1's
    nonNeighbourhood = ones(planeSize(1), planeSize(2), numPlanes);
    
    % Iterate over each point in MarkPoints
    for i = 1:size(MarkPoints, 1)
        % Extract the point coordinates
        x = MarkPoints(i, 1);
        y = MarkPoints(i, 2);
        z = MarkPoints(i, 3);
        
        % Check if the z-coordinate is within the valid range
        if z >= 1 && z <= numPlanes
            % Create a grid of coordinates for the current plane
            [X, Y] = meshgrid(1:planeSize(2), 1:planeSize(1));
            
            % Calculate the squared distance from the point to all positions
            % in the plane
            distances = sqrt((X - x).^2 + (Y - y).^2);
            
            % Set the areas within the specified radius to 0 (neighbourhood)
            nonNeighbourhood(:, :, z) = nonNeighbourhood(:, :, z) & (distances > radius);
        end
    end
    Neighbourhood=abs(nonNeighbourhood-1);
end


