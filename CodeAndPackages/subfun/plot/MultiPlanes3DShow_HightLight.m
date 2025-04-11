function MultiPlanes3DShow_HightLight(Img, cellBoundary, Pos3D, Pos3Dlabel, Zdepth, colorCell, ImgClim,Highlight3D,varargin)
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

    if nargin<9
        PlotParam.EdgeParam=[0.06 0.1 0.06 0.06 0.06 0.06];
        PlotParam.CellCenterWith=1;
        PlotParam.CellBoundaryWidth=0.5;
        PlotParam.PlotCenter=1;
        PlotParam.Alpha=1;
        PlotParam.HighLightColor=[0 1 0];
        PlotParam.HighLightWidth=2;
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

    % If Img is 2D
    if length(Dim) == 2
        % Display the image using imagesc
        imagesc(Img');
        
        % If cell boundaries are provided, plot them
        if ~isempty(cellBoundary)
            plotCellBoundary3D(cellBoundary, [], colorCell, PlotParam.CellBoundaryWidth);
        end
        if ~isempty(Pos3Dlabel)
            labelCellCenter(Pos3D(:,[2 1]), Pos3Dlabel,colorCell);
        end

    % If Img is 3D (multi-plane)
    elseif length(Dim) == 3
        % Loop through each plane and display it in a subplot

         MultiMatrix3DPlotZ(Img,Zdepth,PlotParam.Alpha);
         colormap(gray);
            
            % Adjust the image display contrast using ImgClim
          caxis(ImgClim);
            if ~isempty(cellBoundary)
               plotCellBoundary3D(cellBoundary, Pos3D(:,3), colorCell, PlotParam.CellBoundaryWidth);
            end

            if (~isempty(Pos3D))&PlotParam.PlotCenter
               % plotCellCenter(Pos3D(I,[2 1]), 9, colorCell(I,:),PlotParam.CellCenterWith)
               plotCellCenter3D(Pos3D, 9,colorCell)
            end
            if ~isempty(Pos3Dlabel)
               labelCellCenter(Pos3D, Pos3Dlabel,colorCell);
             end

            plotCellCenter3D(Highlight3D, 12,PlotParam.HighLightColor,PlotParam.HighLightWidth);

            colormap(gray);
            
            % Adjust the image display contrast using ImgClim
            caxis(ImgClim);
            set(gca,'xlim',[0 512],'ylim',[0 512],'zlim',Zdepth([1 end]),'View',[64 24],'zDir','reverse');


            set(gca,'xtick',[],'ytick',[],'ztick',[]);
        end

        % Handle any other cases (if needed)
    end



