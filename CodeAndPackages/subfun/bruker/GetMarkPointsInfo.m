%% 
function [Points,GenInfo,TriggerInfo]=GetMarkPointsInfo(Path)

%Get information from *MarkPoitns.xml files in folder Path

% PathF='E:\LuRecording\11022023SLM\TSeries-11022023-1625-005\';

MPList=dir([Path '*MarkPoints*.xml']);
PointsInfo=struct([]);
for i=1:length(MPList)

    MPstr(i)=xml2struct([MPList(i).folder '\' MPList(i).name]);

    temp=MPstr(i).PVMarkPointSeriesElements;
    temp1=temp.Attributes;
    if ~isfield(temp1,'Category')
       temp1.Category=[];
    end
    if ~isfield(temp1,'Name')
       temp1.Name=[];
    end
    GenInfo(i)=temp1;


    temp2=temp.PVMarkPointElement;
    TriggerInfo(i)=temp2.Attributes;

    tempPoint=temp2.PVGalvoPointElement.Point;
    for iPoint=1:length(tempPoint)
        if length(PointsInfo)==0
           PointsInfo=tempPoint{iPoint}.Attributes;
        else
           PointsInfo(end+1)=tempPoint{iPoint}.Attributes;
            
        end
    end
end

%% Convert struct to Table for easier use
Points=struct2table(PointsInfo);
uniqueRows = unique(Points, 'rows');
Points = convertCellTable(uniqueRows);






