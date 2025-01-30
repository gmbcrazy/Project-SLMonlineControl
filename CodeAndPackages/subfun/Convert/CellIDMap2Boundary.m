% Function to convert cellIDMap to cellBoundary using bwboundaries

function cellBoundary = CellIDMap2Boundary(cellIDMap,cellCandidate)
    
if nargin==1
    % Identify unique non-zero cell IDs
    cellCandidate = setdiff(unique(cellIDMap), 0);
    % Remove NaN values from the candidate list
    cellCandidate(isnan(cellCandidate)) = [];
    cellCandidate=sort(cellCandidate);
end
 
    % Loop through each cell candidate
    for icell = 1:length(cellCandidate)
        % Create a temporary binary map for the current cell
        tempMap = zeros(size(cellIDMap));
        tempMap(cellIDMap == cellCandidate(icell)) = 1;
        
        % Find the boundary of the cell using bwboundaries
        cellBoundary{icell} = bwboundaries(tempMap);
    end

end

