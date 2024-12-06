function plotCellCenter3D(cellCenter, Radius, varargin)
    % Check the number of input arguments
    if nargin < 3
        colorCell = jet(size(cellCenter, 1));
        LineWidth = 1;
    elseif nargin == 3
        colorCell = varargin{1};
        LineWidth = 1;
    elseif nargin == 4
        colorCell = varargin{1};
        LineWidth = varargin{2};
    else
        % Handle other cases if needed
    end

    % If colorCell is a single row, replicate it for each cellIDmark
    if size(colorCell, 1) == 1
        colorCell = repmat(colorCell, size(cellCenter, 1), 1);
    end
    if length(Radius) == 1
        Radius = repmat(Radius, size(cellCenter, 1), 1);
    end

    if size(cellCenter,2)==2
       cellCenter(:,3)=0;
    end

    % Loop through each center point
    for i = 1:size(cellCenter, 1)
        % Extract x, y, and z coordinates of the current center
        % x = cellCenter(i, 2);
        % y = cellCenter(i, 1);
        % z = cellCenter(i, 3);  % Add the z-coordinate

        % Plot a circle around the center with the specified radius, color, and line width
        % rectangle('Position', [x - Radius(i), y - Radius(i), 2 * Radius(i), 2 * Radius(i)], ...
        %     'Curvature', [1, 1], 'EdgeColor', colorCell(i, :), 'LineWidth', LineWidth,'LineStyle',':');
        % 
       theta = linspace(0, 2*pi, 50);
    
    % Coordinates of the circle in 3D space
       x = cellCenter(i, 2) + Radius(i) * cos(theta);
       y = cellCenter(i, 1) + Radius(i) * sin(theta);
       z= zeros(size(x))+cellCenter(i, 3);
       plot3(x, y, z, 'Color', colorCell(i, :), 'LineWidth',LineWidth,'LineStyle',':');


        % Hold on to overlay circles on the same plot
        hold on;
    end

    % Set axis equal for accurate circle visualization
    % axis equal;

    % Add labels or other plot settings if needed
    % xlabel('X-axis');
    % ylabel('Y-axis');
    % zlabel('Z-axis');
    % title('Cell Centers with Circles (3D)');
end
