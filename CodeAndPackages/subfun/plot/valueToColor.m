function colors = valueToColor(values, clim, cmap)
% valueToColor  Map numeric values to RGB colors from a specified colormap.
%
%   colors = valueToColor(values, clim, cmap) returns an N-by-3 matrix where
%   each row is the RGB color corresponding to an element in the input vector
%   "values". The mapping is done according to the color limits "clim" (a 2-element
%   vector [min max]) and the colormap "cmap" (an M-by-3 matrix).
%
%   Example:
%      % Define a colormap, e.g., jet with 256 colors:
%      cmap = jet(256);
%      % Set color limits:
%      clim = [0, 100];
%      % Define some values:
%      values = [0, 25, 50, 75, 100];
%      % Get the corresponding colors:
%      colors = valueToColor(values, clim, cmap);
%
%   See also colormap, jet, parula.

    % Ensure clim has two elements: [minValue maxValue]
    if numel(clim) ~= 2
        error('Clim must be a two-element vector [min max].');
    end
    
    % Clip the values to lie within the limits
    valuesClipped = max(min(values, clim(2)), clim(1));
    
    % Normalize values to the range [0, 1]
    normalized = (valuesClipped - clim(1)) ./ (clim(2) - clim(1));
    
    % Determine the number of colors in the colormap
    nColors = size(cmap, 1);
    
    % Map normalized values to colormap indices.
    % Multiply by (nColors-1) then add 1 so that 0 -> 1 and 1 -> nColors.
    indices = round(normalized * (nColors - 1)) + 1;
    
    % Ensure indices are within the valid range
    indices = min(max(indices, 1), nColors);
    
    % Extract the corresponding RGB values from the colormap
    colors = cmap(indices, :);
end
