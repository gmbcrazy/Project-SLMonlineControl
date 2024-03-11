%% Lu edited from Ana whiskExtract2021 on 10/26/2023
%major update, in case of multi ROI of whisker detection, combine the multi
% ROI into one.


% --- Ana revision 06/16/2022 --- opts.DataRange = 'A2';
% --- Ana revision 05/30/2022 --- xlsread replacement
% --- Ana last revision 09/21/2021 (no change was made since final data
% processing, except for addition of more padding; all correct) ---

function whiskExtract(tempDirWhisk, processedDirWhisk, excelfilename, rowsOfInterest)

disp('Process started')

% --- Directories ---
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
shift = 1;
% param = struct();
% for i = 1:numel(header)
%     if ~isempty(header{i})
%         header{i}=strrep(header{i},' ','');
%         param.(header{i})=text(1+shift:end,i);
%     end
% end

clear i opts header text

if isempty(rowsOfInterest)
    rowsOfInterest = 1:numel(param.Animal);
end

rawFolders = dir([tempDirWhisk '\']);
rawFolders = rawFolders(~strncmpi('.', {rawFolders.name}, 1));


analyzedFolders = dir(processedDirWhisk);
analyzedFolders = analyzedFolders(~strncmpi('.', {analyzedFolders.name}, 1));

clear i

for m = 1:numel(rowsOfInterest)
    
    k = rowsOfInterest(m)-shift;
    
    % currAnalyzedFolder = analyzedFolders(contains({analyzedFolders.name}, param.Animal{k}) == 1 & contains({analyzedFolders.name}, param.Date{k}) == 1 & contains({analyzedFolders.name}, [param.ExpDivision{k} '_Analysed']) == 1);
    currAnalyzedFolder = analyzedFolders(contains({analyzedFolders.name}, param.Animal{k}) == 1 & contains({analyzedFolders.name}, param.Date{k}) == 1);

    
    fileToLoad = dir([currAnalyzedFolder.folder '\' currAnalyzedFolder.name]); % Assumes that there is one file only
    load ([fileToLoad.folder '\' fileToLoad.name])
    
    % currRawFolder = rawFolders(contains({rawFolders.name}, param.Animal{k}) == 1 & contains({rawFolders.name}, param.Date{k}) == 1);
    currRawFolder = rawFolders(contains({rawFolders.name}, param.Date{k}) == 1);
    % currRawFiles = dir([currRawFolder.folder '\' currRawFolder.name '\Camera 1\*.seq']);
    pathGet = [currRawFolder.folder '\' currRawFolder.name];
    path1=dir(pathGet);
    currExp = path1(contains({path1.name}, param.Date{k}) == 1);
    pathGet = [currExp.folder '\' currExp.name '\Camera 1\'];

    currRawFiles = dir([pathGet '\*.seq']);

    whiskTrialsInd = str2num(param.Bseries{k,1});
    
    whiskTrial = NaN(max(totalFrames), numel(whiskTrialsInd));
    whiskStruc = [];
    for t = 1:numel(whiskTrialsInd)
        
        whiskTrialCurr = currRawFiles(endsWith({currRawFiles.name}, strcat('_',string(whiskTrialsInd(t)),'.seq')));
        fieldname = ['whiskFile' num2str(whiskTrialsInd(t))];
       
        % Load frames trial by trial
        disp(['---Loading ' whiskTrialCurr.name '---'])
        if isempty(whiskTrialCurr)
           continue
        end

        [~, framesWhisk] = Norpix2MATLABAna([whiskTrialCurr.folder '\' whiskTrialCurr.name], [], []);
        framesWhisk=double(framesWhisk);   %%Lu added, change unit8 to double        

        tic
        for r = 1:size(posi,1)
            % Crop images using the manually set ROI
            r1=max(round(posi(r,2)),1);
            r2=min(round(posi(r,2))+round(posi(r,4)),size(framesWhisk,1));

            r3=max(round(posi(r,1)),1);
            r4=min(round(posi(r,1))+round(posi(r,3)),size(framesWhisk,2));
           
            framesWhiskTemp = framesWhisk(r1:r2, r3:r4, :);


            ft=diff(framesWhiskTemp,1,3);     %%%This way is faster
            ft2 = sqrt(squeeze(sum(sum(ft.^2))));
            ft2=[nan;ft2(:)];
            whiskTrial(1:length(ft2),t,r) = ft2;
            whiskStruc.(fieldname)(:,:,r) = ft2;



        end
            whiskStruc.(fieldname) = squeeze(nansum(whiskStruc.(fieldname),3));

        toc

    end
    whiskTrial=nansum(whiskTrial,3);
    
    save([fileToLoad.folder '\' fileToLoad.name], 'whiskTrial', 'whiskStruc', '-append')
end