function [cellIDMap,CellPixCount,MedCenter,cellBoundary]=Suite2pCellIDMapFromStat(stat,FovSize)

%% Get cellIDMap from Suite2p processed data; This function is originally one part of
%roiMatchPub.m https://github.com/ransona/ROIMatchPub



% make masks of cells (for longitidinal tracking etc)
cellIDMap = zeros(FovSize);
MedCenter=[];
[m,n]=size(cellIDMap);


for iCell = 1:length(stat)
    % iCell
    cellIDMapTemp = zeros(FovSize);  %In case of overlaped ROI;
    cellID = iCell;
    % roiPix = sub2ind(size(cellMask),stat{cellID}.ypix+int64(Fall.ops.yrange(1))-1,stat{cellID}.xpix+int64(Fall.ops.xrange(1))-1);
    roiPix = sub2ind(size(cellIDMap),stat{cellID}.ypix+1,stat{cellID}.xpix+1);
    cellIDMap(roiPix) = iCell;
    cellIDMapTemp(roiPix) = 1;
    MedCenter(iCell,1)=max(min(round(median(stat{cellID}.ypix)),m),1); 
    MedCenter(iCell,2)=max(min(round(median(stat{cellID}.xpix)),n),1); 
    cellBoundary{iCell} = bwboundaries(cellIDMapTemp);

end

%% 
temp=setdiff(cellIDMap(:),0);
CellPixCount=histcounts(categorical(temp));
