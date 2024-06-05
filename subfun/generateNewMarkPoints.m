function newMarkPoints = generateNewMarkPoints(MarkPoints, radius, numPlanes, numNewPoints, newPointRadius, planeSize, varargin)
    


% First, generate the non-neighbourhood map using the existing function
    
    [Neighbourhood,nonNeighbourhood] = MarkPoint2Neighbourhood(MarkPoints, radius, numPlanes, planeSize);
    % figure;
    % PlotImagePlane(Neighbourhood)
    % figure;
    % PlotImagePlane(nonNeighbourhood)

    if length(size(Neighbourhood))==3
       Neighbourhood=sum(Neighbourhood,3);
       Neighbourhood=repmat(Neighbourhood,1,1,numPlanes);
       nonNeighbourhood=ones(size(Neighbourhood));
       nonNeighbourhood(Neighbourhood>0)=0;
    end

    %% avoid vessel region for nontarget generating
   if nargin==8
       MeanImg=varargin{1};
       vesselTh=varargin{2};
       VesselMask=zeros(size(MeanImg));
    if numPlanes==1
       MeanImg(:,:)=AmpNormalize(MeanImg(:,:),[0.1 99.9]);
    else
       for iplane=1:numPlanes
       MeanImg(:,:,iplane)=AmpNormalize(MeanImg(:,:,iplane),[0.1 99.9]);
       end
    end

       MeanImg=SmoothDecDim3(MeanImg,[1,1]);

<<<<<<< HEAD
       VesselMask(MeanImg<vesselTh)=1;
%        figure;
%        imagesc(VesselMask);
%        PlotImagePlane(VesselMask)

       VesselMask=SmoothDecDim3(VesselMask,[1,1]);
       VesselMask(VesselMask>0.5)=1;
       % VesselMask(VesselMask<=0.5)=0;
       % VesselMask(VesselMask>=0.33)=1;
%        nonNeighbourhood(VesselMask==0)=1;   %%avoid vessel region for nontarget generating
       nonNeighbourhood(VesselMask==1)=1;   %%using vessel region for nontarget generating
      
%        figure;
%        PlotImagePlane(nonNeighbourhood)

   end

%        figure;
%        PlotImagePlane(nonNeighbourhood);
=======
       % VesselMask(MeanImg<vesselTh)=1;  %% Previouslyavoid vessels as non-target region 
       VesselMask(MeanImg>vesselTh)=1;    %% Currently prefer to vessels region as non-target region 05032024
       nonNeighbourhood(VesselMask==1)=0;   %%avoid vessel region for nontarget generating
   end

>>>>>>> 5d9705633a507e90d89ed4694d4d5458dae57738

    count = 0;
    
    % Try to find valid positions for the new points
    while count < numNewPoints
        % Randomly select a plane
        z = randi(numPlanes);
        % Randomly select x and y within the plane dimensions, ensuring margins for the radius
        x = randi([(newPointRadius+1), (planeSize(2)-newPointRadius)]);
        y = randi([(newPointRadius+1), (planeSize(1)-newPointRadius)]);
        
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
