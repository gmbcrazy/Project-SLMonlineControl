% Bug fixed: 2024-3-13, Lu Zhang
% original version messed up data if there are more than 
% two ImageSequence or Zseries in the TimeSeries in PV. The reason is that
% after first ImageSequence/Zseries excuation, there are some redudent data
% left in Buffer, needs to be clean out.
%            


% ========================
% Convert and save raw data on-line
% Lloyd Russell 2016. Optimised on 2017-03-10. Added 'GUI' on 2017-04-11
%

%Lu Zhang, not sure if the following is necessary. I even want to remove
%the header.

% To do:
% ------
% * Fix first run bug - very first tim e running does not connect to PV
% * Header could be improved by containing bitdepth, and number of channels
%
% File header:
% ------------
% First 2 blocks of output file will contain the following information:
% - pixelsPerLine
% - linesPerFrame


% make figure
handles = [];
handles.fig = figure('Name','PrairieLink RawDataStream',...
    'Position',[50 100 320 100], 'MenuBar','none', 'NumberTitle','off',...
    'Color','w');

% add button
handles.GREEN = [0.05 0.85 0.35];
handles.StartButton = uicontrol('Style','Pushbutton',...
    'Position',[10 10 300 40], 'String','Start', 'FontSize',20,...
    'BackgroundColor',handles.GREEN, 'ForegroundColor','w',...
    'Callback',@ClickStart);

% add text
handles.FileNameText = uicontrol('Style','Text',...
    'Position',[10 75 300 15],'BackgroundColor','w', 'String','(Filename)',...
    'FontWeight','Bold', 'Enable','Inactive', 'ButtonDownFcn',@ClickFilename);

handles.ProgressText = uicontrol('Style','Text',...
    'Position',[10 60 300 15],'BackgroundColor','w', 'String','(Progress)');

% store handles in guidata
guidata(handles.fig, handles)


% callback function
function ClickStart(h, e)
    
    % retrieve guidata
    handles = guidata(h);

    % initialise PrairieLink
    pl = actxserver('PrairieLink.Application');
    pl.Connect();
    pl.SendScriptCommands('-DoNotWaitForScans');
    pl.SendScriptCommands('-LimitGSDMABufferSize true 100');
    pl.SendScriptCommands('-StreamRawData true 50');
    pl.SendScriptCommands('-fa 1');  % set frame averaging to 1

    % get acquisition settings
    samplesPerPixel      = pl.SamplesPerPixel();
    pixelsPerLine        = pl.PixelsPerLine();
    linesPerFrame        = pl.LinesPerFrame();
    totalSamplesPerFrame = samplesPerPixel*pixelsPerLine*linesPerFrame;

%     Commented by Lu Zhang 2/23/2024
%     yaml = ReadYaml('settings.yml');
%     flipEvenRows         = yaml.FlipEvenLines;  % toggle whether to flip even or odd lines; 1=even, 0=odd; %
%     Commented by Lu Zhang 2/23/2024

    flipEvenRows         = 1;  % toggle whether to flip even or odd lines; 1=even, 0=odd;   Added by Lu Zhang 2/23/2024

    % get file name
    baseDirectory = pl.GetState('directory', 1);
    tSeriesName   = pl.GetState('directory', 4);
    tSeriesIter   = pl.GetState('fileIteration', 4);
    tSeriesIter   = sprintf('%0.3d', str2double(tSeriesIter));
    filePath      = [baseDirectory, filesep, tSeriesName '-' tSeriesIter];

    % display file name
    completeFileName = [filePath '.bin'];
    decoration = repmat('*', 1, length(completeFileName));
%     disp([decoration; completeFileName; decoration])
    handles.FileNameText.String = completeFileName;
    handles.StartButton.BackgroundColor = [.8 .8 .8];
    
    % open binary file for writing
    fileID = fopen([filePath '.bin'], 'wb');

    % write file header;  Commented by Lu Zhang, 2024-03-13;
    % fwrite(fileID, pixelsPerLine, 'uint16');
    % fwrite(fileID, linesPerFrame, 'uint16');

    % flush buffer
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
    preview = 0;
    if preview
        figure(2);
        Image = imagesc(zeros(linesPerFrame, pixelsPerLine));
        FrameCounter = title('');
        axis off; axis square; axis tight;
    end

    % get data, do conversion, save to file
    icount=1;
    while running   
        % start timer
        tic;    
%         pause(0.1);

        % get raw data stream (timer = ~20ms)
        [samples, numSamplesRead] = pl.ReadRawDataStream(0); 
%         disp([num2str(numSamplesRead) ' of samples detected'])

        % append new data to any remaining old data
        buffer = [buffer samples(1:numSamplesRead)];   

        % extract full frames
        numWholeFramesGrabbed = floor(length(buffer)/totalSamplesPerFrame);
        IncludedInd=1:numWholeFramesGrabbed*totalSamplesPerFrame;
        toProcess = buffer(IncludedInd);

        % clear data from buffer
%         buffer = buffer((numWholeFramesGrabbed*totalSamplesPerFrame)+1:end);
        buffer(IncludedInd)=[];

%         disp([num2str(length(buffer)) ' of samples detected'])
%         Visulize samples recorded in each loop, only use for debugging!
%         Lu Zhang
        figure(3);
        subplot(3,1,1)
        hold on;plot(icount,numSamplesRead,'r.')
        subplot(3,1,2)
        hold on;plot(icount,length(buffer),'r.')

        subplot(3,1,3)
        if exist('numWholeFramesGrabbed')
%         hold on;plot(icount,numWholeFramesGrabbed,'b.')
        icount=icount+1;

           hold on;plot(loopCounter,frameNum,'g.')
           hold on;plot(loopCounter,framesCounter,'r.')


        end
       
        % process the acquired frames (timer = ~5ms)
        if numWholeFramesGrabbed > 0
            for i = 1:numWholeFramesGrabbed
                if started == 0
                    started = 1;
                end

                % get single frame
                frame = toProcess(((i-1)*totalSamplesPerFrame)+1:(i*totalSamplesPerFrame));

                % process the frame (C++ mex code)
                frame = PrairieLink_ProcessFrame(frame, samplesPerPixel, linesPerFrame, pixelsPerLine, flipEvenRows);

                % save processed frames to file
                fwrite(fileID, frame, 'uint16');

                % increment frame counter
                frameNum = frameNum + 1;

                % debugging: preview plot
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
        handles.ProgressText.String = msg;
        drawnow

        % increment counters
        framesCounter = framesCounter + numWholeFramesGrabbed;
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

        % exit loop if finished (if no data collected for previous X loops)
        if started && loopCounter > 100 && sum(allSamplesRead(end-99:end)) == 0
            running = 0;
        elseif started && loopCounter > 10 && sum(allSamplesRead(end-9:end)) == 0   % Keep running but clean buffer during no-data period (such as MarkPoints) but recording not finished yet (if no data collected for previous Y loops)
            buffer=[];
            running = 1;
        else

        end
    end

    % clean up
    pl.Disconnect();
    delete(pl);
    fclose(fileID);
%     fprintf(['\n' 'Finished!' '\n'])
    handles.StartButton.BackgroundColor = handles.GREEN;
end

function ClickFilename(h, e)
    % retrieve guidata
    handles = guidata(h);
    
    % get filename
    CompleteFilePath = handles.FileNameText.String;
   
    % extract path
    [FolderPath,FileName,FileExt] = fileparts(CompleteFilePath);
    
    % open explorer at current file directory
    dos(['explorer ' FolderPath]);
end
