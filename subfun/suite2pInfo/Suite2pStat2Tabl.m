function statTable=Suite2pStat2Tabl(stat)

%% Get cellproperties from Suite2p processed data; Only considered valid ROI;

cellInfoTemp=stat;

for i=1:length(stat)
    if i==1
    cellInfo(i)=cellInfoTemp{i};
    else
    cellInfo=structAdd(cellInfo,cellInfoTemp{i});
    end
end

statTable=struct2table(cellInfo);

