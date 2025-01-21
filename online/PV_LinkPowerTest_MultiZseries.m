
%% Automatic excute multiple MarkPoint.xml files within a specific Folder
%Noted that the TimeSeries in PrairieLink should be MarkPoint current,
%setting -> Multiple Imaging/Zseries Sequence recording of specific fixed
%frame number
% the 1st Imaing/Zseries served as 1st baseline, Not synchronized
%with any MPs, the rest of Zseries is synchronized with MPs, each MPs should be associated with a specific gpl file and a xml files
%: Zseries (50 repetitions) -> Zseries (50 repetitions, syn MP current)-> Zseries (50 repetitions, syn MP current)-> Zseries (50 repetitions, syn MP current)


%% Lu Zhang 2024, developed from PV_LinkExcuteFolder.m, this is for multiple
%MarkPoints and multiple Function Groups defined in the same .gpl file, 
%Excute multiple .xml files, each .xml file define a specific laser power and specific functional group defined in above .gpl file




% callback function
function [XMLTable,FileGenerateInfo]=PV_LinkPowerTest_MultiZseries(XMLparam,PVparam)


     ProcessFolder=XMLparam.ProcessFolder;
     DoRegistration=XMLparam.DoRegistration;
     numPlanes=PVparam.nPlane;

     if DoRegistration    
        if ~isempty(XMLparam.RegRefOps)
           ops=XMLparam.RegRefOps;
           refImg=XMLparam.RegRefImg;  

        end

     end



     SaveDataFolder=[ProcessFolder 'Data\'];
     LogDataFolder=[ProcessFolder '\DataLog\'];

     mkdir(SaveDataFolder);
     mkdir(LogDataFolder);


    %% Initialize PV_Link, looked into the embedded function plIntial.m for details.
    pl = actxserver('PrairieLink.Application');
    pl.Connect();
    [samplesPerPixel,pixelsPerLine, linesPerFrame, totalSamplesPerFrame, flipEvenRows,baseDirectory,tSeriesName,tSeriesIter,tSeriesIterID]=plIntial(pl,SaveDataFolder);


%%     % Determine the source of MarkPoint files; load from structure or directory'
     XMLpattern = 'R(\d+)Laser([\d.]+)GPoint\s?(\d+)';
     ProcessFolder=XMLparam.ProcessFolder;

     Laser=XMLparam.Laser;

     Round=XMLparam.RoundID;
     ProcessFolder=XMLparam.ProcessFolder;

     %%Define the Cells needs to be test for this Tseries excution;
     PointList=XMLparam.PointList;
     if length(Laser)==1&&length(PointList)>1
        Laser=repmat(Laser,1,length(PointList));
     end


     % MarkPointList=dir([ProcessFolder 'R' num2str(Round) 'Laser' num2str(Laser) '*Point*.xml']);
%      MarkPointList=dir([ProcessFolder 'R' num2str(Round) 'Laser*Point*.xml']);
     for iP=1:length(PointList)
         MarkPointList(iP)=dir([ProcessFolder 'R' num2str(Round) 'Laser' num2str(Laser(iP)) 'GPoint' num2str(PointList(iP)) '.xml']);
     end
     

     [roundIDs, AllpointIDs, laserPowers] = XMLPatterExtract(MarkPointList, XMLpattern);
%      MarkPointList=MarkPointList(ismember(AllpointIDs,PointList));
%      [roundIDs, pointIDs, laserPowers] = XMLPatterExtract(MarkPointList, XMLpattern);
     [MarkPointList,roundIDs,pointIDs,laserPowers]=GetXMLFile(MarkPointList,XMLpattern,Round);
     
     PVparam.ZRepetition=31;
     Ziteration=length(PointList)+1;   %%Noted that the 1st Z is not syn with MP.
%      InterMPFrame=Ziteration*PVparam.ZRepetition*PVparam.nPlane;
     InterMPFrame=PVparam.InterMPRepetition*PVparam.nPlane;
     CumInterMPFrame=cumsum(InterMPFrame);
     StartMPFrame=[0 CumInterMPFrame(1:end-1)];
     maxFrame=PVparam.maxFrame;



    %  ShamPossibility=XMLparam.ShamPossibility;
    %  ExXMLList = generateExecutableList(MarkPointList, PVparam.TrialMPSwitch, XMLparam.ShamPossibility);
    % [roundIDs, laserPowers] = FunXMLPatterExtract(ExXMLList, XMLpattern);
     % ExroundIDs=[];
     % ExlaserPowers=[];
%     ExXMLList
    % ExGPLPointList=GPLPointList(roundIDs);


    preview=0;
    BreakYet=0;
    FlushYet=0;
    CheckRedundant=0;
    HavedChecked=0;


     filePath = [baseDirectory, tSeriesName '-' tSeriesIter];
     XMLTable=[repmat(tSeriesIterID,length(pointIDs),1) roundIDs(:) pointIDs(:) laserPowers(:)];
     XMLTable=array2table(XMLTable,'VariableNames',{'FileID','Round','Point','Lasers'});


     completeFileName = filePath;
     binFile=[completeFileName '.bin'];
     logFile=[LogDataFolder  filesep tSeriesName '-' tSeriesIter '.txt'];
     matFile=[baseDirectory 'ExpInfo-' tSeriesIter '.mat']; 

     FileGenerateInfo.FileKey=[tSeriesName '-' tSeriesIter];
     FileGenerateInfo.binFile=binFile;
     FileGenerateInfo.FileID=tSeriesIterID;
     FileGenerateInfo.logFile=logFile;
     FileGenerateInfo.matFile=matFile;
     FileGenerateInfo.tifFolder=[filePath '\'];
     FileGenerateInfo.gplFile={MarkPointList.gplname};
     FileGenerateInfo.xmlFile={MarkPointList.name};
     FileGenerateInfo.checkingTiffBinMatch=0;
     FileGenerateInfo.motionFile=[];
     FileGenerateInfo.motionMed=[];
     FileGenerateInfo.binFileRaw=[];
 


     if DoRegistration
        fileID = fopen(binFile, 'wb'); %Write online motion collection data to a bin file
        shiftsAndCorrFileID = fopen([filePath '_ShiftsAndCorr.bin'],'wb'); %Write online motion frame by frame
        FileGenerateInfo.motionFile=[filePath '_ShiftsAndCorr.bin'];
        % fileIDraw = fopen([filePath 'Raw.bin'], 'wb');%Write raw imaging without motion correction to a bin file
        FileGenerateInfo.binFileRaw=[filePath 'Raw.bin'];
    else
        fileID = fopen(binFile, 'wb');
    end
 

     LogfileID = fopen(logFile,'w'); %Records of key events to a log file.


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
    motionMed = [];
%          loopTimes      = [];
    droppedData    = [];

    ixml=1;
    loadxml=1;
    BreakPointFrame=CumInterMPFrame(ixml);
    SLMChecking=0;

    
    while running   
        % start timer
%               tic;
%           pause(0.02);

            if loadxml==1&&frameNum>StartMPFrame(ixml)+XMLparam.SwitchXMLPostMPFrame*PVparam.nPlane&&ixml<=length(CumInterMPFrame)+0.1   %%when frame number is 10 frames after the previous MP stimuli, update the next xml file.

               if ixml<=length(CumInterMPFrame)-0.9
                  %ExroundIDs(ixml)=roundIDs(ixml);
                  %ExlaserPowers(ixml)=laserPowers(ixml);
                  ixml
                  pl.SendScriptCommands(['-LoadMarkPoints ' MarkPointList(ixml).folder '\' MarkPointList(ixml).gplname ' True'] ); %%Clear existing MarkPoints and load MarkPoints from the gpl file.
                  pause(0.05);
                  pl.SendScriptCommands(['-LoadMarkPoints ' MarkPointList(ixml).folder '\' MarkPointList(ixml).name] );
                  pause(0.01);
%                   [ExGPLPointList(ixml).name ExXMLList{ixml}]
%                   LogMessage(LogfileID,['LoadMarkPoints FunGroup' num2str(ExroundIDs(ixml)) 'with laser' num2str(ExlaserPowers(ixml)) ' at ' num2str(frameNum)]);   
                  LogMessage(LogfileID,['Load' [MarkPointList(ixml).gplname ' ' MarkPointList(ixml).name] ' at ' num2str(frameNum)]);   


               end
                  BreakPointFrame=CumInterMPFrame(ixml);                                     %Update next break point once a MP stimuli was done
                  LogMessage(LogfileID,['BreakPoint Updated Frame ' num2str(BreakPointFrame)]);   
                  loadxml=0;
                  BreakYet=0;
                  SLMChecking=0;
                  framesCounter=frameNum;
                  ZeroSampleSMLCount=0;
                  SampleSML=zeros(4,1);
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
        buffer(IncludedInd)=[];


        % Update the frame, meantime checking if redudant data exist, if it
        % is, remove them
          if numWholeFramesGrabbed > 0
             framesCounter = framesCounter + numWholeFramesGrabbed;
                if BreakYet==0&&framesCounter>=BreakPointFrame
                   LogMessage(LogfileID,['Check Redundant' num2str(frameNum) ' ' num2str(framesCounter)]);
                   BreakYet=1;
                   buffer=[];
                   CheckRedundant=1;
                   SLMChecking=1;
                   loadxml=1;
                   if BreakPointFrame<CumInterMPFrame(end)&&framesCounter==BreakPointFrame
                      LogMessage(LogfileID,'framesCounter just match BreakPointFrame, checking if Tiff and Bin match!');
                      FileGenerateInfo.checkingTiffBinMatch=1;
                   end
%                    disp(['Hi' num2str(frameNum) ' ' num2str(framesCounter)])                   
                end

                for i = 1:numWholeFramesGrabbed
                      if started == 0
                        started = 1;
                      end
                      plane = mod(frameNum,numPlanes)+1;
                % get single frame
                      frame = toProcess(((i-1)*totalSamplesPerFrame)+1:(i*totalSamplesPerFrame));
                      frame = PrairieLink_ProcessFrame(frame, samplesPerPixel, linesPerFrame, pixelsPerLine, flipEvenRows);

                      frameNum = frameNum + 1;
%                       disp([num2str(frameNum) ' ' num2str(framesCounter)])
                      if CheckRedundant==1 && frameNum>BreakPointFrame+0.1       %%Noted that current frameNum should never be higher than the BreakPointFrame (the last frame before a MP stimuli occures). If there is, that frame is redundant. 
                         CheckRedundant=0;
                         frameNum=frameNum-1;                                    %%Redundant frame is found, do not write data to frame, thus substract the frameNumber by 1.
                         LogMessage(LogfileID,['Redundant Detected' num2str(frameNum) ' ' num2str(framesCounter)]);
%                          ixml=ixml+1;
                         break;                                                  %%Redundant frame is found, do not write data to frame, break the loop of updating.
                      elseif CheckRedundant==1 && frameNum<=BreakPointFrame
                         LogMessage(LogfileID,['Keep frame ' num2str(frameNum)]);
                         CheckRedundant=1;
                         if i == numWholeFramesGrabbed
                              CheckRedundant=0;
%                               LogMessage(LogfileID,['No Redundnant detected, and stop checking breakpoint']);
                         end      
                      else 

                      end

                      if DoRegistration
%                       [regFrame,dv,cv] = return_offsets_phasecorr(single(gpuArray(frame)),ops{plane});
                        [regFrame,dv,cv] = return_offsets_phasecorr(single((frame)),ops{plane});
                        motionTemp=sum(abs(dv));
                        motionMed=[motionMed;motionTemp];
                       % save processed frame and correlation values to file
                          fwrite(fileID, uint16(regFrame), 'uint16');
                          fwrite(shiftsAndCorrFileID, [dv cv], 'single');
                          % fwrite(fileIDraw, uint16(frame), 'uint16');

                      else
                          fwrite(fileID, frame, 'uint16');
                      end





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


          %%Check if there is continous iteras with no sample data right
          %%after SLM; if it is, clean buffer;
          if BreakPointFrame<=CumInterMPFrame(end)&&SLMChecking==1
            while SLMChecking==1
                 pause(0.04);
                 [~, SampleTemp] = pl.ReadRawDataStream(0); 
               
                 if SampleTemp==0
                      ZeroSampleSMLCount=ZeroSampleSMLCount+1;
                 else
                      ZeroSampleSMLCount=0;
                 end
             % append new data to any remaining old data
                if ZeroSampleSMLCount>=3
                   LogMessage(LogfileID,[num2str(ZeroSampleSMLCount) ' iterations without samples during SLM period ']);
                   buffer=[];
                   SLMChecking=0;
                end
            end
          end 
          %%Check if there is continous iteras with no sample data right
          %%after SLM; if it is, clean buffer;


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
            subplot(1,3,1)
            hold on;
            plot(loopCounter,numSamplesRead,'r.')
            ylabel('numSamples');
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
%                buffer=[];
%                running = 1;
%                [samples, numSamplesRead] = pl.ReadRawDataStream(0);
%             end

            if started && loopCounter > 600 && sum(allSamplesRead(end-400:end)) == 0 % Keep running but clean buffer during no-data period (such as MarkPoints) but recording not finished yet (if no data collected for previous Y loops)
               LogMessage(LogfileID,[num2str(frameNum) ' frames saved, no more samples detected, terminated.']);
               running=0;
            end
    end

         fclose(fileID);
          if DoRegistration
             fclose(shiftsAndCorrFileID);
             % fclose(fileIDraw);
              motionMed=median(motionMed);
              LogMessage(LogfileID,['Median motion of ' num2str(motionMed) ' pixels (MotionX + MotionY) detected']);
              FileGenerateInfo.motionMed=motionMed;
          end
         fclose(LogfileID);
     save(matFile,'FileGenerateInfo','XMLTable','XMLparam','PVparam');

    %% Update file name for next recording trial



    pl.Disconnect();
    delete(pl);


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





function [samplesPerPixel,pixelsPerLine, linesPerFrame, totalSamplesPerFrame, flipEvenRows,baseDirectory,tSeriesName,tSeriesIter,tSeriesIterID]=plIntial(pl,SaveDataFolder)

    
         pl.SendScriptCommands(['-SetSavePath ' SaveDataFolder]);
         pl.SendScriptCommands('-DoNotWaitForScans');
         pl.SendScriptCommands('-LimitGSDMABufferSize true 100');
         pl.SendScriptCommands('-StreamRawData true 120');
         pl.SendScriptCommands('-fa 1');  % set frame averaging to 1

         samplesPerPixel      = pl.SamplesPerPixel();
         pixelsPerLine        = pl.PixelsPerLine();
         linesPerFrame        = pl.LinesPerFrame();
         totalSamplesPerFrame = samplesPerPixel*pixelsPerLine*linesPerFrame;
        flipEvenRows         = 1;  % toggle whether to flip even or odd lines; 1=even, 0=odd;
    % get file name
        baseDirectory = pl.GetState('directory', 1);
        tSeriesName   = pl.GetState('directory', 4);
        tSeriesIter   = pl.GetState('fileIteration', 4);
%     tempIter=num2str(tSeriesIter)
        tSeriesIter   = sprintf('%0.3d', str2double(tSeriesIter));

  %     pl.SendScriptCommands(['-LoadMarkPoints ' MarkPointGPL.folder '\' MarkPointGPL.name] );

        tSeriesIterID=str2num(tSeriesIter);

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


function [roundIDs, laserPowers] = FunGPLPatterExtract(MarkPointList, GPLpattern)
    % PatterExtract Extracts round IDs, point IDs, and laser powers from a list of filenames.
    % Inputs:
    %   MarkPointList - A structure array containing file information, typically from a dir() call.
    %   XMLpattern - A regular expression pattern designed to extract specific numerical IDs from the filenames.
    % Outputs:
    %   roundIDs - An array containing point identifiers.
    %   laserPowers - An array of laser power values extracted from file names.

    % Initialize arrays to hold extracted data
%     laserPowers = zeros(length(MarkPointList), 1);
    roundIDs = zeros(length(MarkPointList), 1);

    % Iterate through each file in the MarkPointList
    for ixml = 1:length(MarkPointList)
        % Extract the file name from the structure
        fileName = MarkPointList{ixml};
        
        % Use regex to parse out the desired data from the file name
        tokens = regexp(fileName, GPLpattern, 'tokens');
        
        % Check if the regex pattern matched and tokens were found
        if ~isempty(tokens)
            % Convert the captured strings to numbers and store them in the respective arrays
%             laserPowers(ixml) = str2double(tokens{1}{1});
            roundIDs(ixml) = str2double(tokens{1}{1});
        end
    end

    % Round the extracted numeric values to ensure they are integers
    roundIDs = round(roundIDs);
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




