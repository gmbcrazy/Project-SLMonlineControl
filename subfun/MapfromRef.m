%% 
function [Map,CommonCell]=MapfromRef(CrossSessMapping)

%% Cells recorded under the same view across days, try to find the cells consistent recorded from different days.
% 

% Noted that CrossSessMapping is a cell variable. CrossSessMapping{i} is field of 'allSessionMapping' field of
% roiMatchData, where roiMatchData is output variable of function roiMatchPub.m, mapping cells between session i and i+1.

SessionNum=length(CrossSessMapping)+1;
Mapping1=CrossSessMapping;
Map=ForwardMapping(Mapping1);
CommonCell=NormMapWithRef(Map);

disp([num2str(size(CommonCell,1)) ' cells were consistent recorded across ' num2str(size(CommonCell,2)) ' days']);

% if RefSession<=2
%    %ForwardMapping from RefSession; No backwardMapping
%    Map=ForwardMapping(Mapping1);
% elseif RefSession>=SessionNum-1 
%    %BackwardMapping from RefSession; No ForwardMapping
%    Map=BackwardMapping(Mapping1);
% 
% else 
%    %ForwardMapping from RefSession;
%    %BackwardMapping from RefSession; 
%    ForwardMap=Mapping1(RefSession:end);
% 
%    ForwardMap=ForwardMapping(ForwardMap);
%    UpdateMap=[Mapping1(1:RefSession-1) ForwardMap];
%    UpdateMap=BackwardMapping(UpdateMap);
%    Map=ForwardMapping(UpdateMap);
% 
% end

end
%%
function MapBack=BackwardMapping(Mapping2)

Mapping3=Mapping2;
for i=length(Mapping3):-1:2
    [~,i1,i2]=intersect(Mapping3{i}(:,1),Mapping3{i-1}(:,2));
    Mapping3{i}=Mapping3{i}(i1,:);
    Mapping3{i-1}=Mapping3{i-1}(i2,:);
end

Mapping2=Mapping3;
for i=1:length(Mapping2)-1
    [~,i1,i2]=intersect(Mapping2{i}(:,2),Mapping2{i+1}(:,1));
    Mapping2{i}=Mapping2{i}(i1,:);
    Mapping2{i+1}=Mapping2{i+1}(i2,:);
end

MapBack=Mapping2;
end
%%
function MapForward=ForwardMapping(Mapping2)






Mapping2=Mapping2;
for i=1:length(Mapping2)-1
    [~,i1,i2]=intersect(Mapping2{i}(:,2),Mapping2{i+1}(:,1));
    Mapping2{i}=Mapping2{i}(i1,:);
    Mapping2{i+1}=Mapping2{i+1}(i2,:);
end

Mapping3=Mapping2;
for i=length(Mapping3):-1:2
    [~,i1,i2]=intersect(Mapping3{i}(:,1),Mapping3{i-1}(:,2));
    Mapping3{i}=Mapping3{i}(i1,:);
    Mapping3{i-1}=Mapping3{i-1}(i2,:);
end

MapForward=Mapping3;

    

end

%%

function MapNorm=NormMapWithRef(CrossSessMapping)



%%Finished

SessionNum=length(CrossSessMapping)+1;
%%%Noted that CrossSessMapping
MapNorm=ForwardMatch(CrossSessMapping);
% MapNorm=BackwardMatch(CrossSessMapping);

% if RefSession == 1
% 
% end
% 
% if RefSession<=2
%    %ForwardMapping from RefSession; No backwardMapping
%    Map=ForwardMapping(Mapping1);
% elseif RefSession>=SessionNum-1 
%    %BackwardMapping from RefSession; No ForwardMapping
%    Map=BackwardMapping(Mapping1);
% 
% else 
%    %ForwardMapping from RefSession;
%    %BackwardMapping from RefSession; 
%    ForwardMap=Mapping1(RefSession:end);
% 
%    ForwardMap=ForwardMapping(ForwardMap);
%    UpdateMap=[Mapping1(1:RefSession-1) ForwardMap];
%    UpdateMap=BackwardMapping(UpdateMap);
%    Map=ForwardMapping(UpdateMap);
% 
% end

end



%% 
function MapSession=BackwardMatch(Mapping2)

MapSession=Mapping2{end};

for i=length(Mapping2)-1:-1:1
    [~,i2]=ismember(Mapping2{i+1}(:,1),Mapping2{i}(:,2));
    Mapping2{i}=Mapping2{i}(i2,:);
    MapSession=[Mapping2{i}(:,1) MapSession];
end
end



%% 
function MapSession=ForwardMatch(Mapping2)

MapSession=Mapping2{1};

for i=1:length(Mapping2)-1
    [~,i2]=ismember(Mapping2{i}(:,2),Mapping2{i+1}(:,1));
    Mapping2{i+1}=Mapping2{i+1}(i2,:);
    MapSession(:,end+1)=Mapping2{i+1}(:,2);
end
    

end