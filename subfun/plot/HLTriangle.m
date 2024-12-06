function HLTriangle(EndPoint, TrangleRotation, TriangleSize, Color)
    % Check if a single set of parameters is provided for all triangles
    if numel(TriangleSize) == 1
        TriangleSize = repmat(TriangleSize, size(EndPoint, 1), 1);
    end
    if size(Color, 1) == 1
        Color = repmat(Color, size(EndPoint, 1), 1);
    end
    if numel(TrangleRotation) == 1
        TrangleRotation = repmat(TrangleRotation, size(EndPoint, 1), 1);
    end
    hold on;
    % Plot each triangle
    for i = 1:size(EndPoint, 1)
        x = EndPoint(i, 1);
        y = EndPoint(i, 2);
        sizeVal = TriangleSize(i);
        colorVal = Color(i);
        
        % Fixed Cornangles for the corners
        ArrowCorner = pi / 5;
        Cornangles = [ArrowCorner, (pi - ArrowCorner) / 2, (pi - ArrowCorner) / 2];
        
        % Calculate vertices of the triangle with specified rotation
        % vertices = [x, y;...
        %             x + sizeVal * cos(TrangleRotation(i) - pi/2 + Cornangles(1)), y + sizeVal * sin(TrangleRotation(i) - pi/2 + Cornangles(1));...
        %             x + sizeVal * cos(TrangleRotation(i) - pi/2 + Cornangles(2)), y + sizeVal * sin(TrangleRotation(i) - pi/2 + Cornangles(2));...
        %             x + sizeVal * cos(TrangleRotation(i) - pi/2 + Cornangles(3)), y + sizeVal * sin(TrangleRotation(i) - pi/2 + Cornangles(3))];
        vertices = [x, y;...
                    x + sizeVal * cos( - pi/2 + Cornangles(1)), y + sizeVal * sin( - pi/2 + Cornangles(1));...
                    x + sizeVal * cos( - pi/2 + Cornangles(2)), y + sizeVal * sin( - pi/2 + Cornangles(2));...
                    x + sizeVal * cos( - pi/2 + Cornangles(3)), y + sizeVal * sin( - pi/2 + Cornangles(3))];

        % Plot the triangle with no edge color
        fill(vertices(:, 1), vertices(:, 2), colorVal, 'EdgeColor', 'none');
    end
end
