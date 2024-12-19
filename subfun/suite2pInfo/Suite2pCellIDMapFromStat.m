function [cellIDMap, CellPixCount, MedCenter, cellBoundary, varargout] = Suite2pCellIDMapFromStat(stat, FovSize)

%% Get cellIDMap from Suite2p processed data
% This function generates a cell identification map (cellIDMap) using the cell
% statistics (`stat`) provided by Suite2p, a popular pipeline for
% processing calcium imaging data. It assigns unique identifiers to each
% cell and provides information about cell boundaries, cell centers, and pixel counts.
% Originally a part of `roiMatchPub.m` (https://github.com/ransona/ROIMatchPub)

% Initialize an empty cell ID map of the specified field of view (FovSize)
cellIDMap = zeros(FovSize); 
MedCenter = []; % Initialize matrix to hold cell center coordinates
[m, n] = size(cellIDMap); % Dimensions of the field of view (FOV)

% Initialize a 3D cell ID map where each cell is represented separately
cellIDMap3D = repmat(cellIDMap, [1, 1, length(stat)]);

% Loop through each cell in the `stat` structure
for iCell = 1:length(stat)
    % Create a temporary cell ID map for each cell to handle overlapping regions
    cellIDMapTemp = zeros(FovSize); 
    cellID = iCell; % Cell ID corresponds to the loop index

    % Convert the pixel indices of the cell (from `stat`) to linear indices for the FOV
    % `sub2ind` converts the (x, y) coordinates to linear indices
    % Note: `+1` adjusts for MATLAB's 1-based indexing
    roiPix = sub2ind(size(cellIDMap), stat{cellID}.ypix + 1, stat{cellID}.xpix + 1);

    % Assign the cell ID to the corresponding pixels in the map
    cellIDMap(roiPix) = iCell;
    cellIDMapTemp(roiPix) = 1; % Mark the pixels in the temporary cell ID map

    % Calculate the median center for each cell
    MedCenter(iCell, 1) = max(min(round(median(stat{cellID}.ypix)), m), 1); % y-coordinate
    MedCenter(iCell, 2) = max(min(round(median(stat{cellID}.xpix)), n), 1); % x-coordinate

    % Extract the boundaries of the cell using `bwboundaries` for later use
    cellBoundary{iCell} = bwboundaries(cellIDMapTemp);

    % Store the temporary cell ID map in the 3D map for optional output
    cellIDMap3D(:, :, iCell) = cellIDMapTemp;
end

%% Calculate pixel counts for each unique cell ID
% Count the number of pixels associated with each cell ID in `cellIDMap`
temp = setdiff(cellIDMap(:), 0); % Exclude zero entries (background)
CellPixCount = histcounts(categorical(temp)); % Pixel count per cell

% If the user requests the 5th output argument, return the 3D cell ID map
if nargout == 5
   varargout{1} = cellIDMap3D;
end






















%%Original function without comment by ChatGPT



% function [cellIDMap,CellPixCount,MedCenter,cellBoundary,varargout]=Suite2pCellIDMapFromStat(stat,FovSize)
% 
% %% Get cellIDMap from Suite2p processed data; This function is originally one part of
% %roiMatchPub.m https://github.com/ransona/ROIMatchPub
% 
% 
% 
% % make masks of cells (for longitidinal tracking etc)
% cellIDMap = zeros(FovSize);
% MedCenter=[];
% [m,n]=size(cellIDMap);
% 
% cellIDMap3D=repmat(cellIDMap,[1,1,length(stat)]);
% 
% for iCell = 1:length(stat)
%     % iCell
%     cellIDMapTemp = zeros(FovSize);  %In case of overlaped ROI;
%     cellID = iCell;
%     % roiPix = sub2ind(size(cellMask),stat{cellID}.ypix+int64(Fall.ops.yrange(1))-1,stat{cellID}.xpix+int64(Fall.ops.xrange(1))-1);
%     roiPix = sub2ind(size(cellIDMap),stat{cellID}.ypix+1,stat{cellID}.xpix+1);
%     cellIDMap(roiPix) = iCell;
%     cellIDMapTemp(roiPix) = 1;
%     MedCenter(iCell,1)=max(min(round(median(stat{cellID}.ypix)),m),1); 
%     MedCenter(iCell,2)=max(min(round(median(stat{cellID}.xpix)),n),1); 
%     cellBoundary{iCell} = bwboundaries(cellIDMapTemp);
%     cellIDMap3D(:,:,iCell)=cellIDMapTemp;
% end
% 
% %% 
% temp=setdiff(cellIDMap(:),0);
% CellPixCount=histcounts(categorical(temp));
% 
% if nargout==5
%    varargout{1}=cellIDMap3D;
% end
% 

