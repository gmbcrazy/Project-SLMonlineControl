function [TBLmat, TBLInfo] = ExtractExp_FromMat(DataFolder)
    % Function to extract experiment information from .mat files in a given folder
    % 
    % Inputs:
    %   DataFolder - Directory containing the .mat files
    % 
    % Outputs:
    %   TBLmat - Table containing extracted FileID and motionMed data
    %   TBLInfo - Table containing extracted ExpKey and TSeriesENVFile information
    
    DataList = dir(fullfile(DataFolder, 'ExpInfo*.mat'));
    
    % Initialize empty tables
    TBLmat = table();
    TBLInfo = table();
    
    for iFile = 1:length(DataList)
        % DataList(iFile).name
        [TBLmatTemp, TBLInfoTemp] = Mat2ExpInfo(fullfile(DataFolder, DataList(iFile).name));
        TBLmat = [TBLmat; TBLmatTemp];
        TBLInfo = [TBLInfo; TBLInfoTemp];
    end
end

function [MatTBL, TBLInfo] = Mat2ExpInfo(PathFile)
    % Function to load and extract specific fields from a .mat experiment file
    % 
    % Inputs:
    %   PathFile - Full path to the .mat file
    % 
    % Outputs:
    %   MatTBL - Table with extracted FileID and motionMed data
    %   TBLInfo - Table with extracted ExpKey and TSeriesENVFile data
    
    tempLoad = load(PathFile);
    if ~isfield(tempLoad.FileGenerateInfo,'FileID')
        tempLoad.FileGenerateInfo.FileID=str2num(tempLoad.FileGenerateInfo.FileKey(end-2:end));
    end
    % Extract FileID and motionMed
    % tempLoad.FileGenerateInfo.motionMed
    if ~isempty(tempLoad.FileGenerateInfo.motionMed)
      MatTBL = table(tempLoad.FileGenerateInfo.FileID, double(tempLoad.FileGenerateInfo.motionMed), 'VariableNames', {'FileID', 'motionMed'});
    else
      MatTBL = table(tempLoad.FileGenerateInfo.FileID, NaN, 'VariableNames', {'FileID', 'motionMed'});

    end
                   

    % Extract ExpKey and TSeriesENVFile safely
    PVparam = tempLoad.PVparam;
    ExpKey = '';
    TSeriesENVFile = '';

    if isfield(PVparam, 'ExpKey')
        ExpKey = PVparam.ExpKey;
    end

    if isfield(PVparam, 'TSeriesENVFile')
        TSeriesENVFile = PVparam.TSeriesENVFile;
    end 

    TBLInfo = table({ExpKey}, {TSeriesENVFile}, 'VariableNames', {'ExpKey', 'TSeriesENVFile'});
end
