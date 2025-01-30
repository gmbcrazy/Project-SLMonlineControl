function [result, mergedRegion, contourPixels] = check_high_value_center(matrix, threshold_percentage, max_distance, min_region_size, PeakTh, min_merged_region_size, contourMethod)
%CHECK_HIGH_VALUE_CENTER Checks for high-value regions, merges them, and outputs the contour pixels.
%
%   [RESULT, MERGEDREGION, CONTOURPIXELS] = CHECK_HIGH_VALUE_CENTER(MATRIX)
%   returns true if there is at least one connected high-value region that
%   is large enough, close to the center of MATRIX, and the peak value of
%   MATRIX exceeds PeakTh. It also returns the merged indices of all regions
%   that meet the criteria and the contour pixels of the merged region.
%
%   [RESULT, MERGEDREGION, CONTOURPIXELS] = CHECK_HIGH_VALUE_CENTER(MATRIX, THRESHOLD_PERCENTAGE, MAX_DISTANCE, MIN_REGION_SIZE, PeakTh, MIN_MERGED_REGION_SIZE, CONTOURMETHOD)
%   allows customization of the high-value threshold, maximum distance to the center,
%   minimum individual region size, the required minimum peak value, the minimum size of the final merged region,
%   and the method to compute contour pixels.
%
%   Inputs:
%       MATRIX - The heat map matrix to analyze.
%       THRESHOLD_PERCENTAGE - (Optional) Fraction of the peak value to define high-value pixels (default: 0.8).
%       MAX_DISTANCE - (Optional) Maximum allowed distance (in pixels) between the median of a region and the center (default: 0.1 * min(matrix dimensions)).
%       MIN_REGION_SIZE - (Optional) Minimum number of connected pixels for individual regions to be considered (default: 10 pixels).
%       PeakTh - (Optional) Minimum required peak value of the matrix (default: 0).
%       MIN_MERGED_REGION_SIZE - (Optional) Minimum number of pixels for the final merged region (default: 20 pixels).
%       CONTOURMETHOD - (Optional) Method to compute contour pixels: 'perim' or 'boundaries' (default: 'perim').
%
%   Outputs:
%       RESULT - Logical true or false indicating if the merged region meets the criteria.
%       MERGEDREGION - An Nx2 array of [row, col] indices of pixels in the merged region.
%       CONTOURPIXELS - An Mx2 array of [row, col] indices of contour pixels of the merged region.

    % Set default values for optional parameters
    if nargin < 2 || isempty(threshold_percentage)
        threshold_percentage = 0.8; % Default to 80% of the peak value
    end
    if nargin < 3 || isempty(max_distance)
        [m, n] = size(matrix);
        max_distance = 0.1 * min(m, n); % Default to 10% of the smallest matrix dimension
    end
    if nargin < 4 || isempty(min_region_size)
        min_region_size = 10; % Default minimum individual region size
    end
    if nargin < 5 || isempty(PeakTh)
        PeakTh = 0; % Default minimum peak value
    end
    if nargin < 6 || isempty(min_merged_region_size)
        min_merged_region_size = 20; % Default minimum merged region size
    end
    if nargin < 7 || isempty(contourMethod)
        contourMethod = 'perim'; % Default contour method
    end

    % Validate inputs
    if threshold_percentage <= 0 || threshold_percentage > 1
        error('threshold_percentage must be between 0 and 1.');
    end
    if max_distance <= 0
        error('max_distance must be a positive number.');
    end
    if min_region_size <= 0
        error('min_region_size must be a positive integer.');
    end
    if PeakTh < 0
        error('PeakTh must be a non-negative number.');
    end
    if min_merged_region_size <= 0
        error('min_merged_region_size must be a positive integer.');
    end
    if ~ismember(contourMethod, {'perim', 'boundaries'})
        error('contourMethod must be either ''perim'' or ''boundaries''.');
    end

    % Check if the peak value of the matrix exceeds PeakTh
    peakValue = max(matrix(:));
    if peakValue <= PeakTh
        % Peak value is not high enough; return false and empty outputs
        result = false;
        mergedRegion = [];
        contourPixels = [];
        return;
    end

    % Get the size of the matrix
    [m, n] = size(matrix);

    % Calculate center coordinates
    center_row = (m + 1) / 2;
    center_col = (n + 1) / 2;

    % Determine the high-value threshold
    threshold = threshold_percentage * peakValue;

    % Create a binary image of high-value pixels
    binaryImage = matrix > threshold;

    % Identify connected components (regions) in the binary image
    cc = bwconncomp(binaryImage, 8); % 8-connectivity

    % Initialize merged mask
    mergedMask = false(size(matrix));

    % Process each connected region
    for i = 1:cc.NumObjects
        % Get pixel indices of the region
        pixelIdxList = cc.PixelIdxList{i};

        % Check if the region meets the minimum size requirement
        if numel(pixelIdxList) >= min_region_size
            % Convert linear indices to subscripts
            [rows, cols] = ind2sub(size(matrix), pixelIdxList);

            % Compute the median position of the region
            median_row = median(rows);
            median_col = median(cols);

            % Compute the Euclidean distance to the center
            distance = sqrt((median_row - center_row)^2 + (median_col - center_col)^2);

            % Check if the distance is within the specified max_distance
            if distance <= max_distance
                % Update the merged mask
                mergedMask(pixelIdxList) = true;
            end
        end
    end

    % After processing all regions, check if the merged region meets the minimum size
    mergedRegionSize = nnz(mergedMask);
    if mergedRegionSize >= min_merged_region_size
        result = true;
        % Extract the indices of the merged region
        [mergedRows, mergedCols] = find(mergedMask);
        mergedRegion = [mergedRows, mergedCols];

        % Find the contour pixels of the merged region
        switch contourMethod
            case 'perim'
                % Use bwperim to find the perimeter pixels
                contourMask = bwperim(mergedMask);
                [contourRows, contourCols] = find(contourMask);
                contourPixels = [contourRows, contourCols];
            case 'boundaries'
                % Use bwboundaries to get ordered boundary coordinates
                boundaries = bwboundaries(mergedMask);
                if ~isempty(boundaries)
                    contourPixels = boundaries{1}; % First object (since merged regions are combined)
                else
                    contourPixels = [];
                end
        end
    else
        % Merged region does not meet the minimum size requirement
        result = false;
        mergedRegion = [];
        contourPixels = [];
    end
end
