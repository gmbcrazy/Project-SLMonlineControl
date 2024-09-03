function tbl=gpl2Table(gplFile)

s = gpl2struct(gplFile);
PVPoint=s.PVGalvoPointList.PVGalvoPoint;
clear PV
for i=1:length(PVPoint)
    PV(i)=PVPoint{i}.Attributes;
end

temp=struct2table(PV);
tbl = convertTableEntries(temp);