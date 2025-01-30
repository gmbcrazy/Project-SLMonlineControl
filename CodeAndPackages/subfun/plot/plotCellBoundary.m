% Function to plot cell boundaries based on cell IDs

function plotCellBoundary(cellIDMap, cellIDmark, varargin)
    % Check the number of input arguments
    if nargin < 3
        colorCell = jet(length(cellIDmark));
        LineWidth = 3;
    elseif nargin == 3
        colorCell = varargin{1};
        LineWidth = 3;
    elseif nargin == 4
        colorCell = varargin{1};
        LineWidth = varargin{2};
    else
        % Handle other cases if needed
    end

    % If colorCell is a single row, replicate it for each cellIDmark
    if size(colorCell, 1) == 1
        colorCell = repmat(colorCell, length(cellIDmark), 1);
    end

    % Get cell boundaries using the CellIDMap2Boundary function
    cellBoundary = CellIDMap2Boundary(cellIDMap, cellIDmark);
    hold on;
    % Loop through each cell ID
    for i = 1:length(cellIDmark)
        % Loop through each boundary of the current cell ID
        for iB = 1:length(cellBoundary{i})
            % Plot the cell boundary using specified color and LineWidth
            plot(cellBoundary{i}{iB}(:, 2), cellBoundary{i}{iB}(:, 1), 'color',colorCell(i, :), 'LineWidth', LineWidth);
        end
    end
end
