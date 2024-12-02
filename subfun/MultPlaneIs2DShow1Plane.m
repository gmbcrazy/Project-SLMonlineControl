function MultPlaneIs2DShow1Plane(Img, cellBoundary, Pos3D, Pos3Dlabel, Zdepth, PlaneI, colorCell, ImgClim)
    % MultPlaneIs2DShow - Displays 2D images of multi-plane data and overlays
    %                     cell boundaries and positions.
    %
    % Syntax: MultPlaneIs2DShow(Img, cellBoundary, Pos3D, Zdepth, colorCell, ImgClim)
    %
    % Inputs:
    %    Img - Image data, can be 2D or 3D (multi-plane).
    %    cellBoundary - Cell boundary data to overlay on the image.
    %    Pos3D - 3D positions of cells (X, Y, Z).
    %    Zdepth - Z-depth values corresponding to each plane.
    %    colorCell - Colors to use for cell boundaries.
    %    ImgClim - Limits for image display (contrast adjustment).
    
    % Determine the dimensions of the input image and position data
    Dim = size(Img);
    DimPos = size(Pos3D);
    
    % Extract the Z-position (3rd column) from Pos3D if it exists
    if DimPos(2) == 3
        PosZ = Pos3D(:, 3);
    else
        PosZ = [];
    end

    % If Img is 2D
    if length(Dim) == 2
        % Display the image using imagesc
        imagesc(Img');
        
        % If cell boundaries are provided, plot them
        if ~isempty(cellBoundary)
            plotCellBoundary3D(cellBoundary, [], colorCell, 0.5);
        end
        if ~isempty(Pos3Dlabel)
            labelCellCenter(Pos3D(:,[2 1]), Pos3Dlabel,colorCell);
        end

    % If Img is 3D (multi-plane)
    elseif length(Dim) == 3
        % Loop through each plane and display it in a subplot            
            % Display the current plane's image
            imagesc(Img(:,:,PlaneI)');      
            if ~isempty(Pos3D)
            % Find cells that are close to the current Z-depth
               I = find(abs(Pos3D(:,3) - Zdepth(PlaneI)) < 0.1);
            % If cell boundaries are provided, plot them for the current plane
                if ~isempty(cellBoundary)
                   plotCellBoundary3D(cellBoundary(I), [], colorCell(I,:), 0.5);
                end
                if ~isempty(Pos3Dlabel)
                   labelCellCenter(Pos3D(I,[2 1]), Pos3Dlabel(I),colorCell(I,:));
                end
            end
            % Set the colormap to grayscale
            colormap(gray);
            
            % Adjust the image display contrast using ImgClim
            caxis(ImgClim);
            
            
            % Plot the centers of the cells in the current plane
            if ~isempty(Pos3D)
                 plotCellCenter(Pos3D(I,[2 1]), 10, colorCell(I,:),1.5)
            end

    else
        % Handle any other cases (if needed)
    end
end



