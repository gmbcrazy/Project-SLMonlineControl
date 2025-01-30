function [cellIDMap,CellPixCount,MedCenter]=Suite2pCellIDMap(Fall,varargin)

%% Get cellIDMap from Suite2p processed data; This function is originally one part of
%roiMatchPub.m https://github.com/ransona/ROIMatchPub

if nargin==1
   ops=Fall.ops;
else
   ops=varargin{1};
end

cellValid = Fall.iscell(:,1);
% make masks of cells (for longitidinal tracking etc)
cellIDMap = zeros(size(ops.meanImg));
MedCenter=[];
[m,n]=size(cellIDMap);

validCellList = find(cellValid(:,1)==1);
for iCell = 1:length(validCellList)
    % iCell
    cellID = validCellList(iCell);
    % roiPix = sub2ind(size(cellMask),Fall.stat{cellID}.ypix+int64(Fall.ops.yrange(1))-1,Fall.stat{cellID}.xpix+int64(Fall.ops.xrange(1))-1);
    roiPix = sub2ind(size(cellIDMap),Fall.stat{cellID}.ypix+1,Fall.stat{cellID}.xpix+1);
    cellIDMap(roiPix) = iCell;

    MedCenter(iCell,1)=max(min(round(median(Fall.stat{cellID}.ypix)),m),1); 
    MedCenter(iCell,2)=max(min(round(median(Fall.stat{cellID}.xpix)),n),1); 

end

%% 
temp=setdiff(cellIDMap(:),0);
CellPixCount=histcounts(categorical(temp));
