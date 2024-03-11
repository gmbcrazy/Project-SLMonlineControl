%%
% --- Ana revision 06/16/2022 --- opts.DataRange = 'A2';
% --- Ana revision 05/18/2022 --- removed excelSheetNum from input vars
% --- Ana revision 05/12/2022 --- xlsread replacement
% --- Ana revision 09/21/2021 ---

function[] = AlignBatch(processedDir, excelfilename, rowsOfInterest)

% --- Set paths ---
expFiles = dir(processedDir);
expFiles = expFiles(~strncmpi('.', {expFiles.name}, 1));
%expFilesWhisk = dir(processedDirWhisk);%%by Hari
%expFilesWhisk = expFilesWhisk(~strncmpi('.', {expFilesWhisk.name},1));%%by Hari

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

for m = 1:numel(rowsOfInterest)
    
    k = rowsOfInterest(m)-shift;
    
    % --- Get data --
    currExp = expFiles(contains({expFiles.name}, param.Animal{k}) == 1 & contains({expFiles.name}, param.Date{k}) == 1 & contains({expFiles.name}, param.ExpDivision{k}) == 1);
    currExpPlane = currExp(1);
    currExpVars = dir([currExpPlane.folder '\' currExpPlane.name '\*vars.mat']);
%%commented by Hari 49-52
    %currExpWhisk = expFilesWhisk(contains({expFilesWhisk.name}, param.Animal{k}) == 1 & contains({expFilesWhisk.name}, param.Date{k}) == 1 & contains({expFilesWhisk.name}, param.ExpDivision{k}) == 1);
    %currExpPlaneWhisk = currExpWhisk(1);
    %currExpVarsWhisk = dir([currExpPlaneWhisk.folder '\' currExpPlaneWhisk.name '\*.mat']);
    
    clear currExp currExpPlane currExpWhisk currExpPlaneWhisk
    
    % --- Align Ca and Whisk ---
    % Prelast varin = 0 or 1, vRec contains LFP or not, respectively
    % Should create the option to read it directly from excel master
   % Align2021Revised(currExpVars, currExpVarsWhisk, 0);%%original
    Align2021Revised(currExpVars, 0);%%By Hari
end