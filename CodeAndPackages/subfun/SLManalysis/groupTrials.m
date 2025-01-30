%%
% --- Ana revision 06/16/2022 --- added multi-experiment option + autosave
% --- Ana revision 05/30/2022 --- no change
% --- Ana revision 05/18/2022 --- removed excelSheetNum from input vars
% --- Ana revision 05/12/2022 --- xlsread replacement
% --- Ana 11/20/2021 ---

function [groups] = groupTrials(processedDir, excelfilename, rowsOfInterest) %, expType)

% --- Get experiment folders under input directory ---
expFiles = dir(processedDir);
expFiles = expFiles(~strncmpi('.', {expFiles.name}, 1));

% --- Read master excel file with experiment info ---
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
%         header{i} = strrep(header{i},' ','');
%         param.(header{i}) = text(1+shift:end,i);
%     end
% end

clear i opts header text

% --- Cycle through each experiment (i.e. animal, date, exp division) ---
for m = 1:numel(rowsOfInterest)
    
    k = rowsOfInterest(m)-shift;
    
    % --- Load vars ---
    currExp = expFiles(contains({expFiles.name}, param.Animal{k}) == 1 & contains({expFiles.name}, param.Date{k}) == 1 & contains({expFiles.name}, param.ExpDivision{k}) == 1);
    currExpPlane = currExp(1);
    currExpVars = dir([currExpPlane.folder '\' currExpPlane.name '\*vars.mat']);
    
    clear currExp currExpPlane
    
    %load([currExpVars.folder '\' currExpVars.name], 'fTrial', 'caTrialsMask', 'caTrialsMaskLED', 'caTrialsMaskLEDCL')%original
    load([currExpVars.folder '\' currExpVars.name], 'fTrial', 'caTrialsMask')
    % --- Set base save path ---
    [~, baseSave, ~] = fileparts([currExpVars.folder '\' currExpVars.name]);
    kk = strfind(baseSave,'_vars');
    baseSave = baseSave(1:kk);
    
    clear kk behavVars
    
    % --- Define experimental groups ---
   % caTrialsMaskFinal = [caTrialsMask caTrialsMaskLEDCL caTrialsMaskLED];%original
    caTrialsMaskFinal = [caTrialsMask];

    f                 = fields(fTrial);
    groups            = [];
   % ledBoth           = find(sum(caTrialsMaskFinal(:,2:3),2)>0);original
    ledBoth           = find(sum(caTrialsMaskFinal,2)>0);% by Hari

    if isempty(ledBoth)
        groups.spontAll = f(caTrialsMaskFinal(:,1)==0);
        groups.stim     = f(caTrialsMaskFinal(:,1)==1);
        groups.sound    = f(caTrialsMaskFinal(:,1)==2);
    else
        spontGen = find(caTrialsMaskFinal(:,1)==0);
        spontAll = setdiff(spontGen,ledBoth);
        split    = spontAll(spontAll > ledBoth(1));
        stim     = find(caTrialsMaskFinal(:,1)==1);
        sound    = find(caTrialsMaskFinal(:,1)==2);
        led      = find(caTrialsMaskFinal(:,3)==1);
        led2     = find(caTrialsMaskFinal(:,3)==2);
        led3     = find(caTrialsMaskFinal(:,3)==3);
        led4     = find(caTrialsMaskFinal(:,3)==4);
        ledCl    = find(caTrialsMaskFinal(:,2)==1);
        
        % caTrialsMaskFinal(split, 4) = 1;
        
        groups.spont      = f(setdiff(spontAll, split));
        groups.spontAll   = f(spontAll);
        groups.spontBreak = f(split);
        %     if strcmp(expType, 'opto') == 1
        groups.spontLed   = f(setdiff(led,[stim;sound])); % Can add Led2, 3, 4 etc in the future
        groups.spontLedCl = f(setdiff(ledCl,[stim;sound]));
        groups.stim       = f(setdiff(stim,ledBoth));
        groups.stimLed    = f(setdiff(led,[spontGen;sound]));
        groups.stimLedCl  = f(setdiff(ledCl,[spontGen;sound]));
        groups.sound      = f(setdiff(sound,ledBoth));
        groups.soundLed   = f(setdiff(led,[spontGen;stim]));
        groups.soundLedCl = f(setdiff(ledCl,[spontGen;stim]));
        %     end
        %     if strcmp(expType, 'neuromod') == 1
        %         groups.spontDrug1 = f(setdiff(led,[stim;sound]));
        %         groups.spontDrug2 = f(setdiff(ledCl,[stim;sound]));
        %         groups.stim       = f(setdiff(stim,ledBoth));
        %         groups.stimDrug1  = f(setdiff(led,[spontGen;sound]));
        %         groups.stimDrug2  = f(setdiff(ledCl,[spontGen;sound]));
        %         groups.sound      = f(setdiff(sound,ledBoth));
        %         groups.soundDrug1 = f(setdiff(led,[spontGen;stim]));
        %         groups.soundDrug2 = f(setdiff(ledCl,[spontGen;stim]));
        %     end
    end
    % --- Save vars ---
    save([currExpVars.folder '\' baseSave 'groups'], 'groups')
end
end