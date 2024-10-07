function TiffTable = ExpInfoMultiTiffFolder(DataFolder)


TiffFile = dir([DataFolder, '*TSeries*']);
TiffTable=struct2table(TiffFile);
TiffTable=TiffTable(TiffTable.isdir,:);
TiffTable=TiffTable(:,1);
% TiffTable=convertToStrings(TiffTable);

for i=1:size(TiffTable,1)
    tempName=[DataFolder TiffTable.name{i} '\'];
    % [totalRepetitions(i,1), framesAfterStimuli(i,:),laser(i,:),~,~,~,~] = ExpInfoTiffIndiFolder(tempName);
    [totalRepetitions(i,1), framesAfterStimuli(i,:),laser(i,:),~,~,~,~,~,StimID{i,1},StimGPLInd{i,1}] = ExpInfoTiffIndiFolder(tempName);
% Zdepth, ZdepthLaser,cycleID,planeID,files,StimID,StimGPLInd
% totalRepetitions, framesAfterStimuli,StimuliPower,Zdepth, ZdepthLaser,cycleID,planeID,files
end

TiffTable.totalRepetitions=totalRepetitions;
TiffTable.PostSLMframe=framesAfterStimuli;
TiffTable.Laser=laser;
TiffTable.StimID=StimID;
TiffTable.StimGPLInd=StimGPLInd;


% if ~isempty(StimID)
% TiffTable.StimID=StimID{1};
% TiffTable.StimGPLInd=StimGPLInd;
% end

TiffTable = addFileID(TiffTable);



end



function TiffTable = addFileID(TiffTable)
    % Add a FileID column based on the name column in TiffTable
    
    % Initialize an array to store FileIDs
    fileIDList = zeros(height(TiffTable), 1);
    
    % Regular expression pattern to match the '-001', '-002', etc.
    pattern = '-(\d{3})$';
    
    % Loop through each row in the table to extract the FileID
    for i = 1:height(TiffTable)
        % Extract the name
        nameStr = TiffTable.name{i};
        
        % Use regexp to extract the 3-digit number and convert to numeric
        tokens = regexp(nameStr, pattern, 'tokens');
        
        if ~isempty(tokens)
            fileIDList(i) = str2double(tokens{1}{1});
        end
    end
    
    % Add the FileID column to the table
    TiffTable.FileID = fileIDList;
end
