
%% Automatic excute multiple MarkPoint.xml files within a specific Folder
%Noted that the TimeSeries in PrairieLink should be MarkPoint current
%setting -> Imaging or Zseries Sequence recording of specific frames
%Lu Zhang 2024, developed from PV_LinkExcuteFolder.m




% callback function
function varargout=PV_LinkExcuteXMLFunGroup(XMLparam,PVparam,confSet,varargin)

if nargin==4
    PSTHparam=varargin{1};
    if ~isempty(PSTHparam)
        CalPSTH=1;
    else
        CalPSTH=0;
    end
else
    CalPSTH=0;
end


     ProcessFolder=XMLparam.ProcessFolder;

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


%%     % Determine the source of MarkPoint files; load from structure or directory'
     XMLpattern = 'Laser([\d.]+)Group\s?(\d+)';
     % Point=XMLparam.Point;
     % Laser=XMLparam.Laser;
     % Round=XMLparam.RoundID;
     ProcessFolder=XMLparam.ProcessFolder;
%      BreakPointFrame=PVparam.BreakPointFrame;


     GPLPointList=dir([ProcessFolder 'GPLFunGroup.gpl']);
     pl.SendScriptCommands(['-LoadMarkPoints ' GPLPointList(1).folder '\' GPLPointList(1).name] );  %Load the gpl file including all MP as well as all Functional Groups.



%      MarkPointList=dir([ProcessFolder 'Laser' num2str(Laser) '*Group' num2str(Point) '.xml']);
     MarkPointList=dir([ProcessFolder 'Laser*Group*.xml']);

%      x=XMLparam.TrialMPSwitch;
     InterMPFrame=PVparam.InterMPRepetion*PVparam.nPlane;
     CumInterMPFrame=cumsum(InterMPFrame);
     StartMPFrame=[0 CumInterMPFrame(1:end-1)];
     maxFrame=PVparam.maxFrame;

%      XMLparam.ShamPossibility=0.1;
%      XMLparam.SwitchXMLPostMPFrame;
     ShamPossibility=XMLparam.ShamPossibility;
     ExXMLList = generateExecutableList(MarkPointList, PVparam.TrialMPSwitch, XMLparam.ShamPossibility);

     % ExcuteIndex=1:length(MarkPointList);
     % [MarkPointList,roundIDs,pointIDs,laserPowers]=GetXMLFile(MarkPointList,XMLpattern,Round);
    [groupIDs, laserPowers] = FunXMLPatterExtract(ExXMLList, XMLpattern);
     ExgroupIDs=[];
     ExlaserPowers=[];


    preview=1;
    BreakYet=0;
    FlushYet=0;
    CheckRedundant=0;
    HavedChecked=0;
    % pl.SendScriptCommands(['-LoadMarkPoints ' MarkPointList(ixml).folder '\' MarkPointList(ixml).gplname] );
    % pause(0.01);

     filePath = [baseDirectory, filesep, tSeriesName '-' tSeriesIter];
     completeFileName = filePath;
%      completeFileName = [filePath MarkPointList(ixml).name(1:end-4)];

     fileID = fopen([completeFileName '.bin'], 'wb');
     LogfileID = fopen([LogDataFolder  filesep tSeriesName '-' tSeriesIter '.txt'],'w');
     % LogMessage(LogfileID,MarkPointList(ixml).name(1:end-4));

     flushing = 1;
    while flushing
          [samples, numSamplesRead] = pl.ReadRawDataStream(0);
          if numSamplesRead == 0
             flushing = 0;
          end
    end

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
%          loopTimes      = [];
    droppedData    = [];
    ixml=1;
    loadxml=1;
    BreakPointFrame=CumInterMPFrame(ixml);
    PVparam.BreakPointFrame
    while running   
        % start timer
%               tic;

            if loadxml==1&&frameNum>StartMPFrame(ixml)+XMLparam.SwitchXMLPostMPFrame*PVparam.nPlane&&ixml<=length(CumInterMPFrame)   %%when frame number is 10 frames after the previous MP stimuli, update the next xml file.

               if ixml<=length(CumInterMPFrame)-1
                  ExgroupIDs(ixml)=groupIDs(ixml);
                  ExlaserPowers(ixml)=laserPowers(ixml);

                  pl.SendScriptCommands(['-LoadMarkPoints ' ProcessFolder  ExXMLList{ixml}] );
                  LogMessage(LogfileID,['LoadMarkPoints FunGroup' num2str(ExgroupIDs(ixml)) 'with laser' num2str(ExlaserPowers(ixml)) ' at ' num2str(frameNum)]);   
               end
               BreakPointFrame=CumInterMPFrame(ixml);                                     %Update next break point once a MP stimuli was done
               loadxml=0;
               BreakYet=0;
               framesCounter=frameNum;
               ixml=ixml+1;
            end
        % get raw data stream (timer = ~20ms)
        [samples, numSamplesRead] = pl.ReadRawDataStream(0); 

        % append new data to any remaining old data
        buffer = [buffer samples(1:numSamplesRead)];
        IncludedInd=[];
        % extract full frames
        numWholeFramesGrabbed = floor(length(buffer)/totalSamplesPerFrame);
        IncludedInd=1:numWholeFramesGrabbed*totalSamplesPerFrame;
        toProcess = buffer(IncludedInd);

        % clear data from buffer
%         buffer = buffer((numWholeFramesGrabbed*totalSamplesPerFrame)+1:end);
        buffer(IncludedInd)=[];
%         RestFrameN=length(buffer)/numWholeFramesGrabbed/totalSamplesPerFrame;

          if numWholeFramesGrabbed > 0
             framesCounter = framesCounter + numWholeFramesGrabbed;
                if BreakYet==0&&framesCounter>=BreakPointFrame
                   LogMessage(LogfileID,['Check Redundant' num2str(frameNum) ' ' num2str(framesCounter)]);
                   BreakYet=1;
                   buffer=[];
                   CheckRedundant=1;
                   loadxml=1;
%                    disp(['Hi' num2str(frameNum) ' ' num2str(framesCounter)])                   
                end



                for i = 1:numWholeFramesGrabbed
                      if started == 0
                        started = 1;
                      end

                % get single frame
                      frame = toProcess(((i-1)*totalSamplesPerFrame)+1:(i*totalSamplesPerFrame));

                      frame = PrairieLink_ProcessFrame(frame, samplesPerPixel, linesPerFrame, pixelsPerLine, flipEvenRows);

                      frameNum = frameNum + 1;
%                       disp([num2str(frameNum) ' ' num2str(framesCounter)])
                      if CheckRedundant==1 && frameNum>BreakPointFrame+0.1
                         CheckRedundant=0;
                         frameNum=frameNum-1;
                         LogMessage(LogfileID,['Redundant Detected' num2str(frameNum) ' ' num2str(framesCounter)]);
%                          ixml=ixml+1;
                         break;
                      elseif CheckRedundant==1 && frameNum<=BreakPointFrame
                         disp(['Keep frame ' num2str(frameNum)])
                         LogMessage(LogfileID,['Keep frame ' num2str(frameNum)]);
                         CheckRedundant=1;
                         if i == numWholeFramesGrabbed
                              CheckRedundant=0;
                              LogMessage(LogfileID,['No Redundnant detected, and stop checking breakpoint']);
                         end      
                      else

                      end

                      fwrite(fileID, frame, 'uint16');
                      if preview
                         Image.CData = frame';
                         FrameCounter.String = msg;
                         pause(0.00001);
                      end

                               
%                      if preview&& (frameNum == BreakPointFrame+1)
%                         figure(figHandle);
%                         subplot(1,3,1)
%                         imagesc(frame);
%                         axis off; axis square; axis tight;
%                         xlabel('1st frame after breakpoint');
%                      end

               end
          end


             msg = ['Frame: ' num2str(frameNum) ', Loop: ' num2str(loopCounter) ', Sample: ' num2str(totalSamples)];
              loopCounter = loopCounter + 1;
              totalSamples = totalSamples + numSamplesRead;
              allSamplesRead(end+1) = numSamplesRead;

              droppedData(end+1) = pl.DroppedData();
              if droppedData(end)
                 LogMessage(LogfileID,['\n!!! DROPPED DATA AT FRAME ' num2str(frameNum) ' !!!\n']);
              end

        % exit loop if finished if maxFrame frames were collected
          if started && frameNum >= maxFrame
              running = 0;
              LogMessage(LogfileID,[num2str(frameNum) ' frame saved from ' num2str(framesCounter) ' frames detected']);
              if preview==1
              saveas(figHandle,[LogDataFolder  filesep tSeriesName '-' tSeriesIter '.png'],'png');
              end
          end
         if preview==1
            figHandle=figure(2);
            subplot(1,3,2)
            hold on;plot(loopCounter,numWholeFramesGrabbed,'r.')
            ylabel('numWholeFramesGrabbed')
            subplot(1,3,3)
            if exist('numWholeFramesGrabbed')
               hold on;plot(loopCounter,frameNum,'g.')
%                hold on;plot(loopCounter,framesCounter,'r.')
%                plot(loopCounter,BreakPointFrame,'y.')
            end
            set(gca,'ylim',[0 maxFrame+3],'ytick',[0 maxFrame])
            ylabel('Frame # recorded')
         end
%             if BreakYet&&started && loopCounter > 10 && sum(allSamplesRead(end-9:end)) == 0   % Keep running but clean buffer during no-data period (such as MarkPoints) but recording not finished yet (if no data collected for previous Y loops)
%                buffer=[];
%                running = 1;
%                [samples, numSamplesRead] = pl.ReadRawDataStream(0);
%             end

            if started && loopCounter > 150 && sum(allSamplesRead(end-119:end)) == 0   % Keep running but clean buffer during no-data period (such as MarkPoints) but recording not finished yet (if no data collected for previous Y loops)
               running=0;
            end
    end

         fclose(fileID);
         fclose(LogfileID);
         % disp(['xml file of excutionID ' num2str(ExcuteIndex(ixml)) ' in xmlList.csv is completed']);
         % disp(['Pause ' num2str(InterTrialDelay(ixml)) 's for next trial']);
         % % pause(InterTrialDelay(ixml));



    %% Update file name for next recording trial



    pl.Disconnect();
    delete(pl);

    if CalPSTH==1
        PSTHmap(:,:,:,ixml)=PSTHmapCal([completeFileName '.bin'],PSTHparam,confSet);
    end
    if CalPSTH==1
       PSTHmap=squeeze(PSTHmap);
       varargout{1} = [tSeriesIterID XMLparam.Point XMLparam.Laser XMLparam.RoundID PSTHparam.TargetPos];
       varargout{2}=PSTHmap;

       if PSTHparam.Plot==1&&preview==1
           PlaneZ=confSet.ETL+confSet.scan_Z;
         figure
         MultiMatrix3DPlotZ(PSTHmap,PlaneZ,0.9);
         caxis(PSTHparam.Clim);
         Radius=10;
         colormap(PSTHparam.ColorMap);
         set(gca,'xlim',[0 confSet.SLM_Pixels_Y],'ylim',[0 confSet.SLM_Pixels_X],'zlim',PlaneZ([1 end]),'ztick',PlaneZ,'View',[64 24],'zDir','reverse');
         % plotCellCenter3D(SinglePxyz(iPoint,:), Radius, [0 1 0],1.5);
         if isfield(PSTHparam,'TargetPos')
            plotCellCenter3D(PSTHparam.TargetPos, Radius, [0 1 0],1.5);
         end
         papersizePX=[0 0 12 20]
         set(gcf, 'PaperUnits', 'centimeters');
         set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
         saveas(gcf,[LogDataFolder  filesep tSeriesName '-' tSeriesIter 'PSTHmap.png'],'png');
         saveas(gcf,[LogDataFolder  filesep tSeriesName '-' tSeriesIter 'PSTHmap.png'],'png');

       end
    else
        varargout{1}=[];
        varargout{2}=[];

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








function [groupIDs, laserPowers] = FunXMLPatterExtract(MarkPointList, XMLpattern)
    % PatterExtract Extracts round IDs, point IDs, and laser powers from a list of filenames.
    % Inputs:
    %   MarkPointList - A structure array containing file information, typically from a dir() call.
    %   XMLpattern - A regular expression pattern designed to extract specific numerical IDs from the filenames.
    % Outputs:
    %   groupIDs - An array containing point identifiers.
    %   laserPowers - An array of laser power values extracted from file names.

    % Initialize arrays to hold extracted data
    laserPowers = zeros(length(MarkPointList), 1);
    groupIDs = zeros(length(MarkPointList), 1);

    % Iterate through each file in the MarkPointList
    for ixml = 1:length(MarkPointList)
        % Extract the file name from the structure
        fileName = MarkPointList{ixml};
        
        % Use regex to parse out the desired data from the file name
        tokens = regexp(fileName, XMLpattern, 'tokens');
        
        % Check if the regex pattern matched and tokens were found
        if ~isempty(tokens)
            % Convert the captured strings to numbers and store them in the respective arrays
            laserPowers(ixml) = str2double(tokens{1}{1});
            groupIDs(ixml) = str2double(tokens{1}{2});
        end
    end

    % Round the extracted numeric values to ensure they are integers
    groupIDs = round(groupIDs);
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
         if PSTHparam.SmoothSD>0
         PSTHtemp=SmoothDecDim3(PSTHtemp,PSTHparam.SmoothSD);
         end
end




function executableList = generateExecutableList(xmlFilesStruct, x, ShamPossibility)
%GENERATEEXECUTABLELIST Generates a list of XML files to be executed.
%   EXECUTABLELIST = GENERATEEXECUTABLELIST(X) returns a cell array containing
%   the names of XML files to be executed X times. The selection is random with
%   Laser0.5 files having ShamPossibility probability and other files having 1-ShamPossibility probability.

    % Get the list of XML files in the current directory
    % xmlFilesStruct = dir('*.xml');
    xmlFiles = {xmlFilesStruct.name};

    % Separate the XML files into categories
    laser05Files = xmlFiles(contains(xmlFiles, 'Laser0.5'));
    laserNonZeroFiles = xmlFiles(~contains(xmlFiles, 'Laser0.5'));
    
    % Probabilities
    probLaser05 = ShamPossibility; % 10% probability for Laser0.5 files
    probLaserNonZero = 1-ShamPossibility; % 90% probability for Laser1.7 files
    
    % Initialize the output list
    executableList = cell(x, 1);
    
    % Generate the executable list
    for i = 1:x
        P=rand();
        if P < probLaser05
            % Select a random Laser0.5 file
            selectedFile = laser05Files{randi(length(laser05Files))};
        else
            % Select a random Laser1.7 file
            selectedFile = laserNonZeroFiles{randi(length(laserNonZeroFiles))};
        end
        % Add the selected file to the list
        executableList{i} = selectedFile;
    end
end

