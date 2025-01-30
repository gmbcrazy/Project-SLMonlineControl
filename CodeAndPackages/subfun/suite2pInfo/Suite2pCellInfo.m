function cellInfo=Suite2pCellInfo(Fall)

%% Get cellproperties from Suite2p processed data; Only considered valid ROI;

cellValid = find(Fall.iscell(:,1)==1);
cellInfoTemp=Fall.stat(cellValid);

for i=1:length(cellValid)
    if i==1
    cellInfo(i)=cellInfoTemp{i};
    else
    cellInfo=structAdd(cellInfo,cellInfoTemp{i});
    end
end

cellInfo=struct2table(cellInfo);

