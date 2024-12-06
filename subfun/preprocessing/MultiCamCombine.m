%% This function is set to draw the ROI.
% --- Ana revision 06/16/2022 --- opts.DataRange = 'A2';
% --- Ana revision 05/30/2022 --- xlsread replacement
% --- Ana last revision 09/21/2021 (no change was made since final data
% processing, except for addition of more padding; all correct) ---

function MultiCamCombine(tempDirWhisk, processedDirWhisk, excelfilename, rowsOfInterest)

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
    % clear BehavStruc
    k = rowsOfInterest(m)-shift;
    % k = rowsOfInterest(m);
    

    %%Commented by Lu, in this whisk folder, the name of folder include the date information, and the animal name as that in Excel,
    % currExp = expFiles(contains({expFiles.name}, param.Animal{k}) == 1 & contains({expFiles.name}, param.Date{k}) == 1);
 
    %%Commented by Lu, in this whisk folder, the name of folder include the date information, and the animal name as that in Excel,
    currExp = expFiles(contains({expFiles.name}, param.Date{k}) == 1);
    pathGet = [currExp.folder '\' currExp.name];

    disp(['process file ' pathGet]);
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

    currExp = path1(contains({path1.name}, param.Date{k}) == 1 & contains({path1.name}, param.DOB{k}) == 1);
    pathGet = [currExp.folder '\' currExp.name '\'];

    % pathGet = [currExp(1).folder '\'];
    
    
    allFiles = dir([pathGet]);
    allFiles = allFiles(~strncmpi('.', {allFiles.name}, 1));
    BehTrials = allFiles(strncmp('whisk', {allFiles.name}, numel('whisk')));
    BehTrialsInd = str2num(param.Bseries{k,1});
    
    [pathSave] = [processedDirWhisk Date '_' Animal '\'];
    pathSave = processedDirWhisk;

    % mkdir(pathSave)
    
    % Open 10 frames per trial (just to evaluate drifts over time; should
    % preallocate
    framesWhisk = [];
    totalFrames = [];
    for t = 1:numel(BehTrialsInd)
        
        BehTrialCurr = BehTrials(endsWith({BehTrials.name}, strcat('_',string(BehTrialsInd(t)),'_results.mat')));
        CamData(t) = load([BehTrialCurr.folder '\' BehTrialCurr.name]);

    end
    
    load([pathSave '\' SavingName]);
    for t = 1:numel(BehTrialsInd)   
    whiskData=getfield(whiskStruc,['whiskFile' num2str(BehTrialsInd(t))]);
    nS=length(whiskData);
    BehMag=zeros(nS,size(CamData(t).difference_magnitude,2)+1)+nan;
    BehDiffX=zeros(nS,size(CamData(t).difference_magnitude,2))+nan;
    BehDiffY=BehDiffX;
    BehMaxLikelihood=BehDiffX;

    LackFrameNum=nS-max(CamData(t).Camera1_syncFrameNums);
    if LackFrameNum<0
       disp('more behavior samples than whisker samples, check mistakes')
       return
    else
       disp([num2str(LackFrameNum) ' less frames in Behavior detecting than the Whisker or Cam1 video' ])
    end

    BehMag(CamData(t).Camera1_syncFrameNums,1:end-1)=CamData(t).difference_magnitude;
    BehMag(:,end)=whiskData;
    BehDiffX(CamData(t).Camera1_syncFrameNums,:)=CamData(t).difference_x;
    BehDiffY(CamData(t).Camera1_syncFrameNums,:)=CamData(t).difference_y;
    BehMaxLikelihood(CamData(t).Camera1_syncFrameNums,:)=CamData(t).max_likelihoods;

    BehLabel=CamData(t).processed_keypoints;
    BehLabel{end+1}='whisk';
    BehavStruc(t).BehData=BehMag;
    BehavStruc(t).BehLabel=BehLabel;
    BehavStruc(t).BehDiffX=BehDiffX;
    BehavStruc(t).BehDiffY=BehDiffY;
    BehavStruc(t).BehMaxLikelihood=BehMaxLikelihood;
   
    end
    save([pathSave '\' SavingName],'BehavStruc','BehLabel','-append');
    clear BehavStruc CamData
end