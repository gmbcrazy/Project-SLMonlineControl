function newMarkPoints = generateNewMarkPoints(MarkPoints, radius, numPlanes, numNewPoints, newPointRadius, planeSize)
    % First, generate the non-neighbourhood map using the existing function
    [Neighbourhood,nonNeighbourhood] = MarkPoint2Neighbourhood(MarkPoints, radius, numPlanes, planeSize);
    
    if length(size(Neighbourhood))==3
       Neighbourhood=sum(Neighbourhood,3);
       Neighbourhood=repmat(Neighbourhood,1,1,numPlanes);
       nonNeighbourhood=ones(size(Neighbourhood));
       nonNeighbourhood(Neighbourhood>0)=0;
    end

    % Convert nonNeighbourhood to a GPU array
    %nonNeighbourhood = gpuArray(nonNeighbourhood);
    
    % Initialize the new MarkPoints
    % newMarkPoints = zeros(numNewPoints, 3, 'gpuArray');
    count = 0;
    
    % Try to find valid positions for the new points
    while count < numNewPoints
        % Randomly select a plane
        z = randi(numPlanes);
        % Randomly select x and y within the plane dimensions, ensuring margins for the radius
        x = randi([newPointRadius+1, planeSize(2)-newPointRadius]);
        y = randi([newPointRadius+1, planeSize(1)-newPointRadius]);
        
        % Calculate indices for the neighbourhood check
        xRange = max(1, x-newPointRadius):min(planeSize(2), x+newPointRadius);
        yRange = max(1, y-newPointRadius):min(planeSize(1), y+newPointRadius);
        
        % Check if the selected point and its neighbourhood is within the non-neighbourhood area
        if all(all(nonNeighbourhood(yRange, xRange, z)))
            % Create a circular mask
            [X, Y] = meshgrid(xRange, yRange);
            distances = sqrt((X - x).^2 + (Y - y).^2);
            mask = distances <= newPointRadius;
            
            % Update the nonNeighbourhood to exclude the new point's neighbourhood
            nonNeighbourhood(yRange, xRange, z) = nonNeighbourhood(yRange, xRange, z) & ~mask;
            
            % Save the new point
            count = count + 1;
            newMarkPoints(count, :) = gather([x, y, z]); % Gather data back to CPU
        end

    end
    
    % Convert the result back to normal array if necessary
    %newMarkPoints = gather(newMarkPoints);
end
