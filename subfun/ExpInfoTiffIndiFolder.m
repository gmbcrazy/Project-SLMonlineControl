function [totalRepetitions, framesAfterStimuli] = ExpInfoTiffIndiFolder(folderPath)
    % AnalyzeImagingSeries - Analyzes 3D imaging time series and outputs
    %                        total repetitions and frame IDs after stimuli.
    %
    % Syntax: [totalRepetitions, framesAfterStimuli] = AnalyzeImagingSeries(folderPath)
    %
    % Inputs:
    %    folderPath - The path to the folder containing the TIFF and XML files.
    %
    % Outputs:
    %    totalRepetitions - The total number of repetitions (#Tiff / #Planes).
    %    framesAfterStimuli - A cell array of file names for the frames immediately after each MarkPoints.xml.
    
    % Get a list of all files in the folder
    files = dir(fullfile(folderPath, '*.ome.tif'));
    xmlFiles = dir(fullfile(folderPath, '*MarkPoints.xml'));
    
    % Initialize arrays to store cycle numbers and plane numbers
    cycleNumbers = [];
    planeNumbers = [];
    
    % Extract cycle numbers and plane numbers from the TIFF files
    for i = 1:length(files)
        fileName = files(i).name;
        % Extract the cycle number and plane number using regular expressions
        tokens = regexp(fileName, '_Cycle(\d{5})_Ch2_(\d{6})\.ome\.tif', 'tokens');
        if ~isempty(tokens)
         cycleNumbers(end+1) = str2num(tokens{1}{1});
         planeNumbers(end+1) = str2num(tokens{1}{2});
         % disp(['Cycle Number: ', cycleNumber]);
         % disp(['Plane Number: ', planeNumber]);
        % else
        %   disp('No match found');
        end
    end
    cycleNumbers=unique(cycleNumbers);
    % Calculate the number of unique planes
    uniquePlanes = unique(planeNumbers);
    numPlanes = length(uniquePlanes);
    
    % Calculate the total repetitions
    totalRepetitions = length(files) / numPlanes;
    
    % Identify frames after each MarkPoints.xml file
    framesAfterStimuli = [];
    
    for i = 1:length(xmlFiles)
        % Extract the cycle number from the MarkPoints.xml file name
        xmlFileName = xmlFiles(i).name;
        xmlTokens = regexp(xmlFileName, 'Cycle(\d{5})', 'tokens');
        if ~isempty(xmlTokens)
            markCycle = str2double(xmlTokens{1}{1});
            % Find the next cycle number after the MarkPoints cycle
            nextCycle = min(find(cycleNumbers > markCycle));
            % if ~isempty(nextCycle)
            %     % Find the corresponding file for the next cycle and plane 001
            %     nextFile = sprintf('TSeries-%%*%05d_Cy*_%03d.ome.tif', nextCycle, uniquePlanes(1));
            %     nextFileMatch = dir([folderPath, nextFile));
            %     if ~isempty(nextFileMatch)
            %         framesAfterStimuli{end+1} = nextFileMatch.name;
            %     end
            % end
        end
        framesAfterStimuli=[framesAfterStimuli;nextCycle];
    end
    if isempty(framesAfterStimuli)
       framesAfterStimuli=0;
    end
end
