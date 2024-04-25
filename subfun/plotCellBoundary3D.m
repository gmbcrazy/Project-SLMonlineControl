function plotCellBoundary3D(cellBoundary, CellZ,varargin)
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

    nCell=length(cellBoundary);

    % If colorCell is a single row, replicate it for each cellIDmark
    if size(colorCell, 1) == 1
        colorCell = repmat(colorCell, nCell, 1);
    end
    hold on;
    % Loop through each center point
    if isempty(CellZ)
    for i = 1:nCell

        
        for iB = 1:length(cellBoundary{i})
            % Plot the cell boundary using specified color and LineWidth
            plot(cellBoundary{i}{iB}(:, 1), cellBoundary{i}{iB}(:, 2), 'color',colorCell(i, :), 'LineWidth', LineWidth);
        end

    end
    else
    for i = 1:nCell

        
        for iB = 1:length(cellBoundary{i})
            % Plot the cell boundary using specified color and LineWidth
            plot3(cellBoundary{i}{iB}(:, 1), cellBoundary{i}{iB}(:, 2), zeros(size(cellBoundary{i}{iB}(:, 1)))+CellZ(i),'color',colorCell(i, :), 'LineWidth', LineWidth);
        end

    end


    end

end
