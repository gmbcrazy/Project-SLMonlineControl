whiskTrial1%% This function is set to draw the ROI.

% --- Lu revision 11/27/2023 --- adjusting data version of
% develop project, for data without multi camera tracking. 

% --- Ana revision 06/16/2022 --- opts.DataRange = 'A2';
% --- Ana revision 05/30/2022 --- xlsread replacement
% --- Ana last revision 09/21/2021 (no change was made since final data
% processing, except for addition of more padding; all correct) ---

function whiskAdjust_Camera1_2023(tempDirWhisk, processedDirWhisk, excelfilename, nPos,rowsOfInterest)

% --- Directories ---
expFiles = dir(tempDirWhisk);
expFiles = expFiles(~strncmpi('.', {expFiles.name}, 1));

% [~,text,~] = xlsread(excelfilename,excelSheetNum);
opts = detectImportOptions(excelfilename); opts.DataRange = 'A2';
text = readmatrix(excelfilename, opts);

header = opts.VariableNames;
for i = 1:numel(header)
    if ~isempty(header{i})
        header{i} = strrep(header{i},' ','');
        param.(header{i}) = text(1:end,i);
    end
end
% header = text(1,:);

shift = 1;     %%%%What's this param?  This refers to the row in Excel file. Since the 1st row is usually variable names, so not used. The rows of Interest should exclude this row
% shift = 0;     %%%%What's this param?  This refers to the row in Excel file. Since the 1st row is usually variable names, so not used. The rows of Interest should exclude this row


% param = struct();
% for i = 1:numel(header)
%     if ~isempty(header{i})
%         header{i}=strrep(header{i},' ','');
%         param.(header{i})=text(1+shift:end,i);
%     end
% end

clear i opts header text
% rowsOfInterest=[]   %%Lu Added
if isempty(rowsOfInterest)
    rowsOfInterest = 1:numel(param.Animal);
end

for m = 1:numel(rowsOfInterest)      %%%%%%%%%%%Loop for each row, and draw the ROI for whiskers detections.
% for m = 1:length(rowsOfInterest)      %%%%%%%%%%%Loop for each row, and draw the ROI for whiskers detections.
    
    k = rowsOfInterest(m)-shift;
    % k = rowsOfInterest(m);
    

    %%Commented by Lu, in this whisk folder, the name of folder include the date information, and the animal name as that in Excel,
    % currExp = expFiles(contains({expFiles.name}, param.Animal{k}) == 1 & contains({expFiles.name}, param.Date{k}) == 1);
 
    %%Commented by Lu, in this whisk folder, the name of folder include the date information, and the animal name as that in Excel,
    currExp = expFiles(contains({expFiles.name}, param.Date{k}) == 1);
    pathGet = [currExp.folder '\' currExp.name];

    % Date=datevec(currExp.date);
    % Date=num2str(Date(1:3));
    % Date(Date==' ')=[];
    Animal=param.Animal{k};
    Date=param.Date{k};
        if isfield(param,'Date_2') == 1
           SavingName=[Date '_' Animal '_' param.Date_2{k}];
        else
           SavingName=[Date '_' Animal];
        end
        % if isfield(param,'ExpDivision') == 1
        % SavingName = [SavingName '_' param.ExpDivision{k,1}];
        % end

    % Animal=Animal(findstr(Animal,'-')+1:end);
    % SubfolderName=[Date '-mouse' Animal-]

    path1=dir(pathGet);

    % currExp = path1(contains({path1.name}, param.Date{k}) == 1&contains({path1.name}, param.DOB{k}) == 1);
    % currExp = path1(contains({path1.name}, param.Date{k}) == 1);

    % & contains({expFiles.name}, param.DOB{k}) == 1
    pathGet = [currExp.folder '\' currExp.name '\Camera 1\'];
    allFiles = dir([pathGet '\']);
    if isempty(allFiles)
       pathGet = [currExp.folder '\' currExp.name '\'];
       allFiles = dir([pathGet '\']);
    end

    allFiles = allFiles(~strncmpi('.', {allFiles.name}, 1));
    whiskTrials = allFiles(strncmp('whisk', {allFiles.name}, numel('whisk')));

    whiskTrialsInd = [];
    whiskTrialsInd = str2num(param.Bseries{k,1});
    
    [pathSave] = [processedDirWhisk Date '_' Animal '\'];
    pathSave = processedDirWhisk;

    % mkdir(pathSave)
    
    % Open 10 frames per trial (just to evaluate drifts over time; should
    % preallocate
    framesWhisk = [];
    totalFrames = [];
    for t = 1:numel(whiskTrialsInd)
        
        whiskTrialCurr = whiskTrials(endsWith({whiskTrials.name}, strcat('_',string(whiskTrialsInd(t)),'.seq')));
        if ~isempty(whiskTrialCurr)

        [infoTemp, framesWhiskTemp] = Norpix2MATLABAna([whiskTrialCurr.folder '\' whiskTrialCurr.name], 10, []);
        % [infoTemp, framesWhiskTemp] = Norpix2MATLABAna([whiskTrialCurr.folder '\' whiskTrialCurr.name], 500, []);
        % implay(framesWhiskTemp,10);
        CamFrameRate=infoTemp.FrameRate;

        framesWhisk = cat(3, framesWhisk, framesWhiskTemp);
        totalFrames(t) = infoTemp.AllocatedFrames;
        totalFiles(t)=whiskTrialCurr;
        else
        totalFrames(t)=0;
        end
    end
    
    clear framesWhiskTemp
    framesWhisk = double(framesWhisk);
    posi = [];

    if ~isempty(framesWhisk)

    figure; imshow(std(framesWhisk,1,3), [0 20])
    % figure; imshow(mean(framesWhisk,3), [0 20])

    for r = 1:nPos
        h = drawrectangle; % should probably use polygon
        pause
        disp 'Press enter when ready'
        posi(r, :) = h.Position;
    end
    close all
    end

    save([pathSave '\' SavingName], 'framesWhisk', 'totalFrames', 'posi','totalFiles','whiskTrialsInd','CamFrameRate')
    % save([pathSave '\' SavingName], 'CamFrameRate','-append')

    clear totalFiles
end