% Function to add cell labels to based on cell IDs

function labelCellCenter(cellCenter, cellIDmark, varargin)
    % Check the number of input arguments
    if nargin < 3
        colorCell = jet(size(cellCenter,1));
        FontSize = 8;
    elseif nargin == 3
        colorCell = varargin{1};
        FontSize = 8;
    elseif nargin == 4
        colorCell = varargin{1};
        FontSize = varargin{2};
    else
        % Handle other cases if needed
    end

    % If colorCell is a single row, replicate it for each cellIDmark
    if size(colorCell, 1) == 1
        colorCell = repmat(colorCell, size(cellCenter,1), 1);
    end

    if isnumeric(cellIDmark)
       for i=1:length(cellIDmark)
        cellLabel{i}=num2str(cellIDmark(i));
       end
    elseif iscell(cellIDmark)

       cellLabel=cellIDmark;
    else
    end

    % Loop through each center point
    if size(cellCenter,2)==2
        for i = 1:size(cellCenter, 1)
        % Extract x and y coordinates of the current center
        x = cellCenter(i, 2);
        y = cellCenter(i, 1);
        % Plot a circle around the center with the specified radius, color, and line width
        text(x, y, cellLabel{i}, 'Color', colorCell(i, :), 'FontSize', FontSize,'horizontalalignment','center','verticalalignment','base','FontWeight','bold');

        % Hold on to overlay circles on the same plot
        hold on;
        end
    
    elseif size(cellCenter,2)==3
        for i = 1:size(cellCenter, 1)

        x = cellCenter(i, 2);
        y = cellCenter(i, 1);
        z = cellCenter(i, 3);

        % Plot a circle around the center with the specified radius, color, and line width
        text(x, y, z, cellLabel{i}, 'Color', colorCell(i, :), 'FontSize', FontSize,'horizontalalignment','center','verticalalignment','base','FontWeight','bold');

        hold on;
         end
    else


    end
    % Set axis equal for accurate circle visualization
    % axis equal;

    % Add labels or other plot settings if needed
    % xlabel('X-axis');
    % ylabel('Y-axis');
    % title('Cell Centers with Circles');
end


