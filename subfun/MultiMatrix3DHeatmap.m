function MultiMatrix3DHeatmap(Data)
    % Data is a 3D matrix, where each sample is Data(:, :, i)

    % Get the dimensions of the data
    [rows, cols, numMatrices] = size(Data);

    % Create a figure
    figure;

    % Initialize a common set of X and Z coordinates
    X = 1:cols;
    Z = 1:rows;

    % Iterate through the matrices along the Y axis
    for y = 1:numMatrices
        % Extract the current heatmap data
        currentData = squeeze(Data(:, :, y));

        % Plot the heatmap using imagesc
        imagesc(X, Z, currentData);

        % Adjust the colormap and color scaling if needed
        colormap('parula');
        caxis([min(Data(:)), max(Data(:)]);

        % Set labels and title
        xlabel('X-axis');
        ylabel('Z-axis');
        zlabel('Y-axis');
        title(['Matrix ' num2str(y)]);

        % Customize the colorbar if desired
        colorbar;

        % Increment the Y coordinate for the next heatmap
        Z = Z + rows + 1; % Adjust the spacing along the Y axis

        hold on;
    end

    % Adjust the Y-axis limits and make sure it is in the correct order
    ylim([0, (rows + 1) * numMatrices]);
    set(gca, 'YDir', 'reverse');

    % You can further customize the plot as needed.

    hold off;
end
