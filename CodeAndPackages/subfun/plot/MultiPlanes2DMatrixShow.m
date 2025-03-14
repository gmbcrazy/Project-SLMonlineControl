function H=MultiPlanes2DMatrixShow(Img, cellBoundary, Pos3D, Pos3Dlabel, Zdepth, colorCell, ImgClim,varargin)
    % MultiPlanes2DShow - Displays 2D images of multi-plane data and overlays
    %                     cell boundaries and positions.
    %
    % Syntax: MultiPlanes2DShow(Img, cellBoundary, Pos3D, Zdepth, colorCell, ImgClim)
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

    if nargin<8
        PlotParam.RowPlot=1;
        PlotParam.RowColNum=1;
        PlotParam.RowColID=1;
        PlotParam.EdgeParam=[0.06 0.1 0.06 0.06 0.06 0.06];
        PlotParam.CellCenterWith=1;
        PlotParam.CellBoundaryWidth=0.5;
        PlotParam.PlotCenter=1;
    else
        PlotParam=varargin{1};
    end


    if size(colorCell,1)==1
       colorCell=repmat(colorCell,DimPos(1),1);

    end
    % Extract the Z-position (3rd column) from Pos3D if it exists
    if DimPos(2) == 3
        PosZ = Pos3D(:, 3);
    else
        PosZ = [];
    end

for i=1:size(Img,4)

        % Loop through each plane and display it in a subplot
        for iplane = 1:size(Img,3)
            H(i,iplane)=subplotLU(size(Img,4), size(Img,3), i, iplane, PlotParam.EdgeParam);
            % subplotLU(3, 4, iplane, 3,  PlotParam.EdgeParam);
            % Display the current plane's image
            imagesc(Img(:,:,iplane,i)');
            
            if ~isempty(Pos3D)
            % Find cells that are close to the current Z-depth
               I = find(abs(Pos3D(:,3) - Zdepth(iplane)) < 0.1);
            
            % If cell boundaries are provided, plot them for the current plane
                if ~isempty(cellBoundary)
                   plotCellBoundary3D(cellBoundary(I), [], colorCell(I,:), PlotParam.CellBoundaryWidth);
                end
                if ~isempty(Pos3Dlabel)
                   labelCellCenter(Pos3D(I,[2 1]), Pos3Dlabel(I),colorCell(I,:));
                end
            end
            colormap(gray);
            
            % Adjust the image display contrast using ImgClim
            caxis(ImgClim(i,:));
            
            
            % Plot the centers of the cells in the current plane
            if (~isempty(Pos3D))&PlotParam.PlotCenter
                 plotCellCenter(Pos3D(I,[2 1]), 9, colorCell(I,:),PlotParam.CellCenterWith)
            end
            set(gca,'xtick',[],'ytick',[]);
        end
end


end



