function [PixShiftFile, Files] = PixShiftLoad(WorkFolder)
    motCorrBinFile = dir(fullfile(WorkFolder, '*_ShiftsAndCorr.bin'));
    if isempty(motCorrBinFile)
        error('No binary files found in the specified folder.');
    end
    
    Files = {motCorrBinFile.name};
    PixShiftFile = zeros(1, numel(motCorrBinFile));
    
    for iFile = 1:numel(motCorrBinFile)
        tempbin = fopen(fullfile(WorkFolder, motCorrBinFile(iFile).name), 'r');
        PixShift = [];
        
        while ~feof(tempbin)  % Check for end of file
            data = fread(tempbin, 2, 'single');
            if numel(data) < 2  % Break if less than expected data is read
                break;
            end
            PixShift(end+1) = sum(data);
            temp = fread(tempbin, 1, 'single');
        end
        
        fclose(tempbin);
        
        if ~isempty(PixShift)
            PixShiftFile(iFile) = median(PixShift);
        else
            PixShiftFile(iFile) = NaN; % Handle empty file case
        end
    end
end
