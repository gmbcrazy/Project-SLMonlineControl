function ColorOut = Value2Color(Values, ColorMapC, Clim)
%VALUE2COLOR maps numeric values to colors based on given colormap and limits.
%
% Inputs:
% Values     - Numeric array of values to be mapped to colors.
% ColorMapC  - Colormap array (N x 3), each row is an RGB color.
% Clim       - Two-element vector specifying the color limits [min max].
%
% Output:
% ColorOut   - RGB color array corresponding to input values.

% Clip values outside the color limits
C1 = Values(:);
C1(C1 > Clim(2)) = Clim(2);
C1(C1 < Clim(1)) = Clim(1);

% Number of colors in colormap
ClimN = size(ColorMapC, 1);

% Create edges for discretization
E1 = linspace(Clim(1), Clim(2), ClimN + 1);

% Map values to color indices
[ColorC1I, ~] = discretize(C1, E1);

% Assign corresponding colors
ColorOut = ColorMapC(ColorC1I, :);

end