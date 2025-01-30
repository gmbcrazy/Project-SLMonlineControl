function HLArrow(EndPoint, Angle, Length, Color)
    % Check if Angle, Length, and Color are scalar or arrays
    if isscalar(Angle)
        Angle = repmat(Angle, size(EndPoint, 1), 1);
    end
    if isscalar(Length)
        Length = repmat(Length, size(EndPoint, 1), 1);
    end

    % If Color is a scalar, repeat it for each arrow
    if isscalar(Color)
        Color = repmat(Color, size(EndPoint, 1), 1);
    end

    % If Color is a 3-digit RGB value, convert it to a matrix
    if numel(Color) == 3 && isnumeric(Color)
        Color = repmat(Color, size(EndPoint, 1), 1);
    end
   
    hold on;
    % Plot each arrow
    for i = 1:size(EndPoint, 1)
        % Calculate arrow components
        u = cosd(Angle(i)) * Length(i);
        v = sind(Angle(i)) * Length(i);

        % Plot the arrow with specified color
        quiver(EndPoint(i, 1) - u, EndPoint(i, 2) - v, u, v, 'Color', Color(i, :), 'LineWidth', 2);
    end

end
