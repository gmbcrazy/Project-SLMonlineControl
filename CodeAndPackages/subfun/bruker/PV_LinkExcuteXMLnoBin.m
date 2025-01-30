
%% Automatic excute multiple MarkPoint.xml files within a specific Folder
%Noted that the TimeSeries in PrairieLink should be MarkPoint current
%setting -> Imaging or Zseries Sequence recording of specific frames
%Lu Zhang 2024, developed from PV_LinkExcuteFolder.m




% callback function
function PV_LinkExcuteXMLnoBin(XMLparam,PVparam,confSet,varargin)

if nargin==4
    PSTHparam=varargin{1};
    if ~isempty(PSTHparam)
        CalPSTH=1;
    else
        CalPSTH=0;
    end
else
    CalPSTH=0
end



%%     % Determine the source of MarkPoint files; load from structure or directory'
     XMLpattern = 'R(\d+)Laser([\d.]+)GPoint\s?(\d+)';
     Point=XMLparam.Point;
     Laser=XMLparam.Laser;
     Round=XMLparam.RoundID;
     ProcessFolder=XMLparam.ProcessFolder;

     maxFrame=PVparam.maxFrame;
     BreakPointFrame=PVparam.BreakPointFrame;


     MarkPointList=dir([ProcessFolder 'R' num2str(Round) 'Laser' num2str(Laser) '*Point' num2str(Point) '.xml']);
     % ExcuteIndex=1:length(MarkPointList);
     [MarkPointList,roundIDs,pointIDs,laserPowers]=GetXMLFile(MarkPointList,XMLpattern,Round);
     SaveDataFolder=[ProcessFolder 'Data\'];
     LogDataFolder=[ProcessFolder '\DataLog\'];


    mkdir(SaveDataFolder);
    mkdir(LogDataFolder);



    % InterTrialDelay = random('uniform',RandomDelayInterval(1),RandomDelayInterval(2),1,length(MarkPointList));

%     pl.SendScriptCommands('-LoadMarkPoints F:\LuSLMOnlineTest\02152024\WriteFileUpdate3\GPL.gpl');
%     pl.SendScriptCommands('-LoadMarkPoints F:\LuSLMOnlineTest\02152024\WriteFileUpdate3\Point5.xml');



    pl = actxserver('PrairieLink.Application');
    pl.Connect();
    pl.SendScriptCommands(['-SetSavePath ' SaveDataFolder]);
    pl.SendScriptCommands('-DoNotWaitForScans');
    pl.SendScriptCommands('-LimitGSDMABufferSize true 100');
%     pl.SendScriptCommands('-StreamRawData true 50');   
    pl.SendScriptCommands('-StreamRawData true 50');
    pl.SendScriptCommands('-fa 1');  % set frame averaging to 1

    samplesPerPixel      = pl.SamplesPerPixel();
    pixelsPerLine        = pl.PixelsPerLine();
    linesPerFrame        = pl.LinesPerFrame();
    totalSamplesPerFrame = samplesPerPixel*pixelsPerLine*linesPerFrame;
    % yaml = ReadYaml('settings.yml');
    % flipEvenRows         = yaml.FlipEvenLines;  % toggle whether to flip even or odd lines; 1=even, 0=odd;
    flipEvenRows         = 1;  % toggle whether to flip even or odd lines; 1=even, 0=odd;
    % get file name
    baseDirectory = pl.GetState('directory', 1);
    tSeriesName   = pl.GetState('directory', 4);
    tSeriesIter   = pl.GetState('fileIteration', 4);
%     tempIter=num2str(tSeriesIter)
    tSeriesIter   = sprintf('%0.3d', str2double(tSeriesIter));

%     pl.SendScriptCommands(['-LoadMarkPoints ' MarkPointGPL.folder '\' MarkPointGPL.name] );

    tSeriesIterID=str2num(tSeriesIter);
%     tSeriesIter   = sprintf('%0.3d', str2double(tSeriesIterID))



    for ixml=1:length(MarkPointList) %%Randomized the MarkPoint Stimulation Order
        % ixml=StimList(jxml);
    % flipEvenRows         = 0;  % toggle whether to flip even or odd lines; 1=even, 0=odd;
        close all
         BreakYet=0;
         FlushYet=0;
         CheckRedundant=0;
         HavedChecked=0;
         pl.SendScriptCommands(['-LoadMarkPoints ' MarkPointList(ixml).folder '\' MarkPointList(ixml).gplname] );
         pause(0.2);
         pl.SendScriptCommands(['-LoadMarkPoints ' MarkPointList(ixml).folder '\' MarkPointList(ixml).name] );
         pause(0.2);

         filePath = [baseDirectory, filesep, tSeriesName '-' tSeriesIter];
         completeFileName = [filePath MarkPointList(ixml).name(1:end-4)];

         
         LogfileID = fopen([LogDataFolder  filesep tSeriesName '-' tSeriesIter '.txt'],'w');
         LogMessage(LogfileID,MarkPointList(ixml).name(1:end-4));

         pause(0.1);

         flushing = 1;
         while flushing
            [samples, numSamplesRead] = pl.ReadRawDataStream(0);
            if numSamplesRead == 0
                flushing = 0;
            end
         end

    % start the current t-series
         pl.SendScriptCommands('-TSeries');
    % initialise state variables, buffer, and counters/records
         running        = 1;
         started        = 1;
         loopCounter    = 1;
         totalSamples   = 0;
         framesCounter  = 0;
         frameNum       = 0;
         buffer         = [];
         allSamplesRead = [];
         msg            = [];
         loopTimes      = [];
         droppedData    = [];



        while running   
        % start timer
              tic;    

        % get raw data stream (timer = ~20ms)
              [samples, numSamplesRead] = pl.ReadRawDataStream(0); 
        % append new data to any remaining old data

              loopCounter = loopCounter + 1;
              totalSamples = totalSamples + numSamplesRead;
              allSamplesRead(end+1) = numSamplesRead;
              loopTimes(end+1) = toc;

        % test for dropped data
              droppedData(end+1) = pl.DroppedData();
              if droppedData(end)
                 LogMessage(LogfileID,['\n!!! DROPPED DATA AT FRAME ' num2str(frameNum) ' !!!\n']);
              end
%             if loopCounter > 150
%             sum(allSamplesRead(end-119:end))
%             end
%             allSamplesRead(end)
         
            if started && loopCounter > 150 && sum(allSamplesRead(end-119:end)) == 0   % Keep running but clean buffer during no-data period (such as MarkPoints) but recording not finished yet (if no data collected for previous Y loops)
               running=0;
            end
         end


         fclose(LogfileID);



    %% Update file name for next recording trial
         tSeriesIterID=tSeriesIterID+1;
         tSeriesIter   = sprintf('%0.3d', tSeriesIterID);

    end

    pl.Disconnect();
    delete(pl);

end








function [MarkPointList,roundIDs,pointIDs,laserPowers]=GetXMLFile(MarkPointList,XMLpattern,Round)


    [roundIDs,pointIDs,laserPowers]=XMLPatterExtract(MarkPointList,XMLpattern);
    NeedXMLIndex=ismember(roundIDs,Round);
    roundIDs=roundIDs(NeedXMLIndex);
    pointIDs=pointIDs(NeedXMLIndex);
    laserPowers=laserPowers(NeedXMLIndex);
    MarkPointList=MarkPointList(NeedXMLIndex);

    %%Find the corresponding gpl files defining MarkPoint based on each xml
    %%files
    if ~isfield(MarkPointList,'gplname')
        for ixml=1:length(MarkPointList)
            MarkPointGPL(ixml)=dir([MarkPointList(ixml).folder '\R' num2str(roundIDs(ixml)) 'GPoint' num2str(pointIDs(ixml)) '.gpl']);
            MarkPointList(ixml).gplname=MarkPointGPL(ixml).name;
        end
    end


 %% Random execution of all XML files to prevent consecutive executions of the same point
    % Initialize the shuffled execution order
%     executionOrder = zeros(length(MarkPointList), 1);
% 
% % Generate a randomized execution order with the constraint
%     while ~isempty(find(executionOrder == 0, 1))
%     % Shuffle the indices
%     shuffledIndices = randperm(length(MarkPointList));
% 
%     % Check for consecutive PointIDs
%     validShuffle = true;
%         for i = 2:length(shuffledIndices)
%              if pointIDs(shuffledIndices(i)) == pointIDs(shuffledIndices(i-1))
%                 validShuffle = false;
%                  break;
%             end
%          end
% 
%     % If the shuffle is valid, use it as the execution order
%          if validShuffle
%             executionOrder = shuffledIndices;
%         end
%     end


%%    MarkPointGPL = MarkPointGPL(executionOrder);
    % MarkPointList = MarkPointList(executionOrder);
    % MarkPointGPL=MarkPointGPL(executionOrder);
    % roundIDs=roundIDs(executionOrder);
    % pointIDs=pointIDs(executionOrder);
    % laserPowers=laserPowers(executionOrder);

%%

% %     RestFiletable=struct2table(MarkPointList);
% % %     GPLtable=struct2table(MarkPointGPL);
% % %     RestFiletable.gplname=GPLtable.name;
% %     RestFiletable.excutionID=[1:length(MarkPointList)]';
% %     writetable(RestFiletable,[MarkPointList(1).folder '\xmlListCurrent.csv'])



end




function [roundIDs, pointIDs, laserPowers] = XMLPatterExtract(MarkPointList, XMLpattern)
    % PatterExtract Extracts round IDs, point IDs, and laser powers from a list of filenames.
    % Inputs:
    %   MarkPointList - A structure array containing file information, typically from a dir() call.
    %   XMLpattern - A regular expression pattern designed to extract specific numerical IDs from the filenames.
    % Outputs:
    %   roundIDs - An array containing numerical round IDs extracted from the file names.
    %   pointIDs - An array containing point identifiers.
    %   laserPowers - An array of laser power values extracted from file names.

    % Initialize arrays to hold extracted data
    roundIDs = zeros(length(MarkPointList), 1);
    laserPowers = zeros(length(MarkPointList), 1);
    pointIDs = zeros(length(MarkPointList), 1);

    % Iterate through each file in the MarkPointList
    for ixml = 1:length(MarkPointList)
        % Extract the file name from the structure
        fileName = MarkPointList(ixml).name;
        
        % Use regex to parse out the desired data from the file name
        tokens = regexp(fileName, XMLpattern, 'tokens');
        
        % Check if the regex pattern matched and tokens were found
        if ~isempty(tokens)
            % Convert the captured strings to numbers and store them in the respective arrays
            roundIDs(ixml) = str2double(tokens{1}{1});
            laserPowers(ixml) = str2double(tokens{1}{2});
            pointIDs(ixml) = str2double(tokens{1}{3});
        end
    end

    % Round the extracted numeric values to ensure they are integers
    roundIDs = round(roundIDs);
    pointIDs = round(pointIDs);
end


function LogMessage(LogFileID,Message)
         disp(Message);
         fprintf(LogFileID, '%s\r\n', Message);
%          fprintf(LogFileID, '\n');
end






function PSTHtemp=PSTHmapCal(BinFile,PSTHparam,confSet)


         PreData = Suite2pSingleChBin2Frame(BinFile, confSet.SLM_Pixels_Y, confSet.SLM_Pixels_X, length(confSet.ETL), PSTHparam.PreInd);
         PostData= Suite2pSingleChBin2Frame(BinFile, confSet.SLM_Pixels_Y, confSet.SLM_Pixels_X, length(confSet.ETL), PSTHparam.PostInd);
         PSTHtemp=squeeze(mean(PostData,3)-mean(PreData,3));
         PSTHtemp=SmoothDecDim3(PSTHtemp,PSTHparam.SmoothSD);

end

