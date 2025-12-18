function [tbl,StimID,StimPower,PowerWeight]=MPxml2Table(MPxmlFile)

s = xml2struct(MPxmlFile);
PVPointInfo=s.PVMarkPointSeriesElements.PVMarkPointElement.PVGalvoPointElement;
clear PV
GroupInfo=PVPointInfo.Attributes;
PVPoint=PVPointInfo.Point;
StimID=GroupInfo.Points;
StimPower=s.PVMarkPointSeriesElements.PVMarkPointElement.Attributes.UncagingLaserPower;

PowerWeight=ones(length(PVPoint),1);

if isfield(s.PVMarkPointSeriesElements.PVMarkPointElement.Attributes,'CustomLaserPercent')
   PowerWeight=s.PVMarkPointSeriesElements.PVMarkPointElement.Attributes.CustomLaserPercent;
end

if iscell(PVPoint)
    for i=1:length(PVPoint)
        PV(i)=PVPoint{i}.Attributes;
    end
else
    for i=1:length(PVPoint)
         PV(i)=PVPoint(i).Attributes;
    end
end

temp=struct2table(PV);
tbl = convertTableEntries(temp);
tbl;