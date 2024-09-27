function [totalRepetitions, framesAfterStimuli,StimuliPower,Zdepth, ZdepthLaser,cycleID,planeID,files] = ExpInfoTiffIndiFolder(folderPath)
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
    MPxmlFiles = dir(fullfile(folderPath, '*MarkPoints.xml'));
    ImgxmlFile = dir(fullfile(folderPath, '*TSeries*.xml'));

    for ixml=1:length(ImgxmlFile)

    if isempty(findstr(ImgxmlFile(ixml).name,'Voltage'))
       yaml=xml2yaml([ImgxmlFile(ixml).folder '\' ImgxmlFile(ixml).name]);
       break
    end
    end
    tifInfo=imfinfo([files(1).folder '\' files(1).name]);
    tempnumPlanes=length(tifInfo);
    % Initialize arrays to store cycle numbers and plane numbers
    cycleNumbers = [];
    planeNumbers = [];

    clear Zdepth
    Zdepth(:,1)=yaml.Zdepth_ETL;

    I1=find(abs(yaml.Zdepth_ETL-yaml.scan_Z(2))<0.01);
    Zdepth(I1,2)=yaml.scan_Z(1);
    
    I2=setdiff(1:length(yaml.Zdepth_ETL),I1);
    for j=1:length(I2)
        Zdepth(I2(j),2)=Zdepth(I2(j),1)-Zdepth(I1,1)+Zdepth(I1,2);
    end

    Zdepth=Zdepth(:,2);
    ZdepthLaser=yaml.Zdepth_laserPower;

    % Extract cycle numbers and plane numbers from the TIFF files
    for i = 1:length(files)
        fileName = files(i).name;
        % Extract the cycle number and plane number using regular expressions
        tokens = regexp(fileName, '_Cycle(\d{5})_Ch2_(\d{6})\.ome\.tif', 'tokens');
        if ~isempty(tokens)
         cycleNumbers(end+1) = str2num(tokens{1}{1});
         planeNumbers(end+1) = max(str2num(tokens{1}{2}),tempnumPlanes);
         % disp(['Cycle Number: ', cycleNumber]);
         % disp(['Plane Number: ', planeNumber]);
        % else
        %   disp('No match found');
        end
    end
    cycleID=cycleNumbers;
    planeID=planeNumbers;
    cycleNumbers=unique(cycleNumbers);
    % Calculate the number of unique planes
    uniquePlanes = unique(planeNumbers);
    numPlanes = length(uniquePlanes);
    
    % Calculate the total repetitions
    if tempnumPlanes>1
       totalRepetitions = length(files);

    else
        totalRepetitions = length(files) / numPlanes;
    end
    % Identify frames after each MarkPoints.xml file
    framesAfterStimuli = [];
    StimuliPower=[];
    for i = 1:length(MPxmlFiles)
        % Extract the cycle number from the MarkPoints.xml file name
        xmlFileName = MPxmlFiles(i).name;
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
        StimuliPower=[StimuliPower;MPxml2yaml([MPxmlFiles(i).folder '\' MPxmlFiles(i).name])];
        
    end
    if isempty(framesAfterStimuli)
       framesAfterStimuli=0;
    end
    if isempty(StimuliPower)
       StimuliPower=NaN;
    end


end
