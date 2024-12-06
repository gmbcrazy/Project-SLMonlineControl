function [IndDayCell,CellIDNest,CellDayStrU]=MultiDay2Slice(Data,Mapping)

%% Input:
% Data is structure variables including suite2p data from across days, Data(1) is simply loading day 1 of suite2p data folder
%Mapping is the 1st output of function MapfromRef.m

%% Output:
% CellDayStrU is a universal cell list for all cells recorded across days. Not necessary each cell were consistent recorded across multi-days. 
% IndDayCell is cell variable, IndDayCell{i} index of cells in day i in the universal cell list of CellDayStrU 
% CellIDNest is cell variable, CellIDNest{i} is index of cells in day i in the iscell list in the Data(i).iscell


%% include all cells, creat cellID and DayID matrix.
CellAll=[];
CellDay=[];
% CellInDay=[];
for iDay=1:length(Data)
    tempC=[1:sum(Data(iDay).iscell(:,1))]';
    for iCell=1:length(tempC)
        % CellAll{end+1,1}=['D' num2str(iDay) 'N' num2str(iCell)];
        CellAll(end+1,iDay)=iCell;
        CellDay(end+1,iDay)=iDay;
    end
end
CellAll0=CellAll;
CellDay0=CellDay;


for iDay=1:length(Data)-1
    i1=Mapping{iDay}(:,1);
    i2=Mapping{iDay}(:,2);

    % % for iDayBack=length(Mapping):-1:iDay
    % % [~,CommonCell]=MapfromRef(Mapping(iDay:end))
    % % end

    iCurrent=find(CellDay(:,iDay)>0.5);
    iNext=find(CellDay(:,iDay+1)>0.5);
    
    for j=1:length(i1)

        %Find cells repeated recorded in previous day and current day.
        %Adding Day ID and Cell ID of current day to their previous day
        j1=find(CellAll(iCurrent,iDay)==i1(j));
        CellDay(iCurrent(j1),iDay+1)=iDay+1;
        CellAll(iCurrent(j1),iDay+1)=i2(j);

        
        %Adding Day ID and Cell ID of previous day to current day 
        j2=find(CellAll(iNext,iDay+1)==i2(j));
        % CellDay(iNext(j2),iDay)=iDay;
        % CellDay(iNext(j2),iDay)=iDay;

        CellDay(iNext(j2),1:iDay)=repmat(CellDay(iCurrent(j1(1)),1:iDay),length(j2),1);
        CellAll(iNext(j2),1:iDay)=repmat(CellAll(iCurrent(j1(1)),1:iDay),length(j2),1);

        CellDay(iNext(j2),1:iDay)=CellDay(iCurrent(j1(1)),1:iDay);
        CellAll(iNext(j2),1:iDay)=CellAll(iCurrent(j1(1)),1:iDay);



    end

    % CellDay(iNext(i2),iDay)=iDay;
    % CellAll(iNext(i2),iDay)=i1;
  
end


%% Using day and cell information to creat a unique universal cell list for all cells recorded 
clear CellDayStr
for iCell=1:size(CellDay,1)
    DayInd=find(CellDay(iCell,:)>0);
    tempStr='';
    for iDay=1:length(DayInd)
        tempStr=[tempStr 'D' num2str(DayInd(iDay)) 'N' num2str(CellAll(iCell,DayInd(iDay)))];
    end
    CellDayStr{iCell,1}=tempStr;
    % CellDayStr{iCell,1}(CellDayStr{iCell,1}=='0')='';
end
length(unique(CellDayStr));

[CellDayStrU,i1]=unique(CellDayStr);
CellDayU=CellDay(i1,:);
CellAllU=CellAll(i1,:);

DayRecord=[];
for iCell=1:length(CellDayStrU)
    i1=findstr(CellDayStrU{iCell},'D');
    i2=findstr(CellDayStrU{iCell},'N');
    DayRecord(iCell)=str2num(CellDayStrU{iCell}(i1(1)+1:i2(1)-1));
end
[DayRecord,i1]=sort(DayRecord);
CellDayU=CellDayU(i1,:);
CellAllU=CellAllU(i1,:);
CellDayStrU=CellDayStrU(i1);

clear numCell
for iDay=1:length(Data)
    numCell(iDay,1)=sum(Data(iDay).iscell(:,1));
    numCell(iDay,2)=sum(unique(CellAll(:,iDay))>0);
end

clear IndDayCell CellIDNest
for iDay=1:length(Data)
    IndDayCell{iDay}=find(CellAllU(:,iDay)>0);
    CellIDNest{iDay}=CellAllU(IndDayCell{iDay},iDay);
end