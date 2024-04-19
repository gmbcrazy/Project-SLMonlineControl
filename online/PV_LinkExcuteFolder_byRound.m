
%% Automatic excute multiple MarkPoint.xml files within a specific Folder
%Noted that the TimeSeries in PrairieLink should be MarkPoint current
%setting -> Imaging or Zseries Sequence recording of specific frames
%Lu Zhang 2024, developed from PV_LinkExcuteFolder.m




% callback function
function PV_LinkExcuteFolder_byRound(ProcessFolder,RandomDelayInterval,maxFrame,BreakPointFrame,Round,ExcuteIndex)

%%     % Determine the source of MarkPoint files; load from structure or directory'
    XMLpattern = 'R(\d+)Laser([\d.]+)GPoint\s?(\d+)';

    if ~isempty(findstr(ProcessFolder,'Current.csv'))
       csv=dir(ProcessFolder);
       SaveDataFolder=[csv.folder '\Data\'];

       MarkPointList=table2struct(readtable(ProcessFolder));
       MarkPointList=MarkPointList(ExcuteIndex);
      [roundIDs,pointIDs,laserPowers]=XMLPatterExtract(MarkPointList,XMLpattern);
    else
       MarkPointList=dir([ProcessFolder '*Point*.xml']);
       ExcuteIndex=1:length(MarkPointList);
       [MarkPointList,roundIDs,pointIDs,laserPowers]=GetXMLFile(MarkPointList,XMLpattern,Round);
       SaveDataFolder=[ProcessFolder 'Data\'];

       % MarkPointGPL=dir([ProcessFolder 'SingleP\*.gpl']);
    end

    mkdir(SaveDataFolder);



    InterTrialDelay = random('uniform',RandomDelayInterval(1),RandomDelayInterval(2),1,length(MarkPointList));

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
         pl.SendScriptCommands(['-LoadMarkPoints ' MarkPointList(ixml).folder '\' MarkPointList(ixml).gplname] );
         pause(0.2);
         pl.SendScriptCommands(['-LoadMarkPoints ' MarkPointList(ixml).folder '\' MarkPointList(ixml).name] );
         pause(0.1);

         filePath = [baseDirectory, filesep, tSeriesName '-' tSeriesIter];
         completeFileName = [filePath MarkPointList(ixml).name(1:end-4)];


    % open binary file for writing
         fileID = fopen([completeFileName '.bin'], 'wb');
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
         started        = 0;
         loopCounter    = 1;
         totalSamples   = 0;
         framesCounter  = 0;
         frameNum       = 0;
         buffer         = [];
         allSamplesRead = [];
         msg            = [];
         loopTimes      = [];
         droppedData    = [];

    % preview image window (only use for debugging!)
         preview = 1;
         if preview
            figure(2);
            subplot(1,3,1)
            Image = imagesc(zeros(linesPerFrame, pixelsPerLine));
            FrameCounter = title('');
            axis off; axis square; axis tight;
         end

    % get data, do conversion, save to file
        while running   
        % start timer
              tic;    

        % get raw data stream (timer = ~20ms)
              [samples, numSamplesRead] = pl.ReadRawDataStream(0); 

        % append new data to any remaining old data
               buffer = [buffer samples(1:numSamplesRead)];

        % extract full frames
        numWholeFramesGrabbed = floor(length(buffer)/totalSamplesPerFrame);
        IncludedInd=1:numWholeFramesGrabbed*totalSamplesPerFrame;
        toProcess = buffer(IncludedInd);

        % clear data from buffer
%         buffer = buffer((numWholeFramesGrabbed*totalSamplesPerFrame)+1:end);
        buffer(IncludedInd)=[];
        RestFrameN=length(buffer)/numWholeFramesGrabbed/totalSamplesPerFrame;



        % process the acquired frames (timer = ~5ms)

          if numWholeFramesGrabbed > 0
             framesCounter = framesCounter + numWholeFramesGrabbed;

                for i = 1:numWholeFramesGrabbed
                      if started == 0
                        started = 1;
                      end

                % get single frame
                      frame = toProcess(((i-1)*totalSamplesPerFrame)+1:(i*totalSamplesPerFrame));

                % process the frame (C++ mex code)
                      frame = PrairieLink_ProcessFrame(frame, samplesPerPixel, linesPerFrame, pixelsPerLine, flipEvenRows);

                % save processed frames to file
                % increment frame counter
                % BreakPointFrame is number of frames in total before a
                % MarkPoint stimuli is applied, PV has redudant data/frame
                % after this; So when this happen, clean the buffer, and do
                % not write this redudant data into .bin file
                      frameNum = frameNum + 1;
                      if BreakYet==0&&frameNum>BreakPointFrame&&framesCounter>BreakPointFrame
                         frameNum=frameNum-1;
                         BreakYet=1;
                         buffer=[];
                         [samples, numSamplesRead] = pl.ReadRawDataStream(0);                         
                         break;
                      end

                      fwrite(fileID, frame, 'uint16');
                      if preview
                         Image.CData = frame';
                         FrameCounter.String = msg;
                         pause(0.00001);
                     end
               end
          end

        % display progress
%         fprintf(repmat('\b', 1, length(msg)));  % delete previous 'message'
             msg = ['Frame: ' num2str(frameNum) ', Loop: ' num2str(loopCounter) ', Sample: ' num2str(totalSamples)];
%         fprintf(msg);
%         handles.ProgressText.String = msg;
%         drawnow

        % increment counters
              % framesCounter = framesCounter + numWholeFramesGrabbed;
              loopCounter = loopCounter + 1;
              totalSamples = totalSamples + numSamplesRead;
              allSamplesRead(end+1) = numSamplesRead;
              loopTimes(end+1) = toc;

        % test for dropped data
              droppedData(end+1) = pl.DroppedData();
              if droppedData(end)
                 fprintf(2, ['\n!!! DROPPED DATA AT FRAME ' num2str(framesCounter) ' !!!\n'])
                 fprintf(msg)
              end



%         Visulize samples recorded in each loop, only use for debugging!Lu Zhang


        % exit loop if finished if maxFrame frames were collected
            if started && frameNum >= maxFrame
                 running = 0;
            end
         if preview
            figure(2);
            subplot(1,3,2)
            hold on;plot(loopCounter,numWholeFramesGrabbed,'r.')
            ylabel('numWholeFramesGrabbed')
            subplot(1,3,3)
            if exist('numWholeFramesGrabbed')
               hold on;plot(loopCounter,frameNum,'g.')
               hold on;plot(loopCounter,framesCounter,'r.')
%                plot(loopCounter,BreakPointFrame,'y.')
            end
            set(gca,'ylim',[0 maxFrame+3],'ytick',[0 maxFrame])
            ylabel('Frame # recorded')
         end
%             if BreakYet&&started && loopCounter > 10 && sum(allSamplesRead(end-9:end)) == 0   % Keep running but clean buffer during no-data period (such as MarkPoints) but recording not finished yet (if no data collected for previous Y loops)
% %                  buffer=[];
%                running = 1;
%                [samples, numSamplesRead] = pl.ReadRawDataStream(0);
%             end

            if started && loopCounter > 150 && sum(allSamplesRead(end-119:end)) == 0   % Keep running but clean buffer during no-data period (such as MarkPoints) but recording not finished yet (if no data collected for previous Y loops)
               running=0;
            end
         end

    % clean up
         fclose(fileID);
         disp(['xml file of excutionID ' num2str(ExcuteIndex(ixml)) ' in xmlList.csv is completed']);
         disp(['Pause ' num2str(InterTrialDelay(ixml)) 's for next trial']);
         pause(InterTrialDelay(ixml));

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
    executionOrder = zeros(length(MarkPointList), 1);

% Generate a randomized execution order with the constraint
    while ~isempty(find(executionOrder == 0, 1))
    % Shuffle the indices
    shuffledIndices = randperm(length(MarkPointList));
    
    % Check for consecutive PointIDs
    validShuffle = true;
        for i = 2:length(shuffledIndices)
             if pointIDs(shuffledIndices(i)) == pointIDs(shuffledIndices(i-1))
                validShuffle = false;
                 break;
            end
         end
    
    % If the shuffle is valid, use it as the execution order
         if validShuffle
            executionOrder = shuffledIndices;
        end
    end


%%    MarkPointGPL = MarkPointGPL(executionOrder);
    MarkPointList = MarkPointList(executionOrder);
    MarkPointGPL=MarkPointGPL(executionOrder);
    roundIDs=roundIDs(executionOrder);
    pointIDs=pointIDs(executionOrder);
    laserPowers=laserPowers(executionOrder);

%%

    RestFiletable=struct2table(MarkPointList);
%     GPLtable=struct2table(MarkPointGPL);
%     RestFiletable.gplname=GPLtable.name;
    RestFiletable.excutionID=[1:length(MarkPointList)]';
    writetable(RestFiletable,[MarkPointList(1).folder '\xmlListCurrent.csv'])



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

