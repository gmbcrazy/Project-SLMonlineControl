clear all
% Load the .mat file that contains the FOVUpdate structure
load('D:\Project1-LocalProcessing\Step1\07-Oct-2025FOV.mat')

% Loop through each field of view (FOV) in the FOVUpdate structure
for i = 1:length(FOVUpdate)

    % Extract the data folder path for the current FOV
    tempPath = FOVUpdate(i).DataFolder;

    % ---- Identify specific trial indices ----
    % Find the index for spontaneous activity:
    % This corresponds to entries where both Point and Group fields are NaN
    sponI = find((isnan(FOVUpdate(i).subT1.Point) & isnan(FOVUpdate(i).subT1.Group)) > 0);
    sponI = sponI(1);  % Take the first spontaneous index

    % Find indices for PowerTest trials:
    % These are trials where Point is not NaN (indicating a defined test)
    PowerTestI = find(~isnan(FOVUpdate(i).subT1.Point) > 0);
    PowerTestI = [PowerTestI(1) PowerTestI(end)];  % Take the first and last test trials

    % Find indices for Group trials:
    % These are trials where Group is not NaN (indicating group-defined tests)
    GroupI = find(~isnan(FOVUpdate(i).subT1.Group) > 0);
    GroupI = [GroupI(1) GroupI(end)];  % Take the first and last group trials

    % Combine all indices of interest into one array
    CheckI = [sponI PowerTestI GroupI];

    % Retrieve the corresponding file keys for the selected trials
    FileCheck = FOVUpdate(i).subT1.FileKey(CheckI);

    % ---- Compare scan parameters across selected XML files ----
    for j = 1:length(FileCheck)
        % Construct the full path to the XML file for each selected trial
        FileTempCheck = [tempPath FileCheck{j} '\' FileCheck{j} '.xml'];

        % Convert XML file to YAML format for easy data access
        yaml(j) = xml2yaml(FileTempCheck);

        % Compute absolute differences in scan amplitude (X and Y)
        % relative to the first file in the set
        DiffX(j) = sum(abs(yaml(j).ScanAmp_X - yaml(1).ScanAmp_X));
        DiffY(j) = sum(abs(yaml(j).ScanAmp_Y - yaml(1).ScanAmp_Y));
    end

    % ---- Check for significant differences ----
    % If there is any deviation in scan amplitude greater than 0.001 units,
    % display the FOV path and the measured differences
    if max(DiffY) > 0.001 || max(DiffX) > 0.001
        suite2pFOVPath{i}  % Display the corresponding FOV path
        DiffY               % Display Y-axis differences
        DiffX               % Display X-axis differences
    end

    % Clear temporary variables before moving to the next FOV
    clear yaml DiffY DiffX
end
