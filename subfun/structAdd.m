
function structData=structAdd(struct1,struct2)

%%%%Add struct2 to struct1
copyFields = intersect(fieldnames(struct2),fieldnames(struct1));

structData=struct1;
AddN=length(struct2);
for iAdd=1:AddN
    L=length(structData);
    for ifield=1:length(copyFields)
        Temp=getfield(struct2,{iAdd},copyFields{ifield});
        structData=setfield(structData,{L+1},copyFields{ifield},Temp);
    end
end


