function HighlightPlane(xLimits,yLimits,ZHighLight,varargin)


    if nargin < 4
        Color = jet(length(ZHighLight));
        LineWidth = 1;
    elseif nargin == 4
        Color = varargin{1};
        LineWidth = 1;
    elseif nargin == 5
        Color = varargin{1};
        LineWidth = varargin{2};
    else
        % Handle other cases if needed
    end

    if size(Color, 1) == 1
        Color = repmat(Color, length(ZHighLight), 1);
    end

hold on;
for iPlot=1:length(ZHighLight)
     plot3DRectangle(xLimits, yLimits,ZHighLight(iPlot),Color(iPlot,:),LineWidth);
end

% rectangle('Position', [x - Radius(i), y - Radius(i), 2 * Radius(i), 2 * Radius(i)], ...
 %     'Curvature', [1, 1], 'EdgeColor', Color(i, :), 'LineWidth', LineWidth,'LineStyle',':');

end




function plot3DRectangle(xLimits, yLimits,z,zcolor,LineWidth)
    % Check if the input limits are valid
    if numel(xLimits) ~= 2 || numel(yLimits) ~= 2
        error('Input limits must be a 2-element array for both x and y.');
    end
    
    % Set the common Z coordinate
    % z = 0;

    % Define the corners of the rectangle
    corners = [xLimits(1), yLimits(1), z;
               xLimits(2), yLimits(1), z;
               xLimits(2), yLimits(2), z;
               xLimits(1), yLimits(2), z;
               xLimits(1), yLimits(1), z]; % Closing the rectangle

    % Plot the rectangle
    plot3(corners(:,1), corners(:,2), corners(:,3), 'color',zcolor,'LineStyle','-', 'LineWidth', 2);

end
        