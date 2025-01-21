% Bug fixed: 2024-3-13, Lu Zhang
% original version messed up data if there are more than 
% two ImageSequence or Zseries in the TimeSeries in PV. The reason is that
% after first ImageSequence/Zseries excuation, there are some redudent data
% left in Buffer, needs to be clean out.
%            


% ========================
% Convert and save raw data on-line
% Lloyd Russell 2016. Optimised on 2017-03-10. Added 'GUI' on 2017-04-11
% Henry Dalgleish 2018 implemented online motion correction (using Suite2P functions)

% Lu Zhang, remove the header, fixed bugs above from version of Lloyd Russell PrairieLink_RawDataStream
% and add Henry Dalgleish motion correction version from PrairieLink_RawDataStreamReg

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
    'Position',[50 100 400 240],... % Adjust height for added space
    'MenuBar','none', 'NumberTitle','off', 'Color','w');
% add button
handles.GREEN = [0.05 0.85 0.35];
handles.StartButton = uicontrol('Style','Pushbutton',...
    'Position',[50 20 300 40], 'String','Start', 'FontSize',20,...
    'BackgroundColor',handles.GREEN, 'ForegroundColor','w',...
    'Callback',@ClickStart);

% add text
handles.DoRegistration = uicontrol('Style','Checkbox',...
    'Position',[20 190 150 20],'BackgroundColor','w', 'String','Do registration');

handles.LoadRefImgButton = uicontrol('Style','Pushbutton',...
    'Position',[200 190 180 20],'BackgroundColor','w', 'String','Load reference image',...
    'Callback',@LoadReferenceImage);

handles.RefImgText = uicontrol('Style','Text',...
    'Position',[20 140 360 40],... % Adjust height for multiline text
    'BackgroundColor','w', 'String',' ',...
    'FontAngle','Italic', 'HorizontalAlignment','left');

handles.FileNameText = uicontrol('Style','Text',...
    'Position',[20 120 360 40],'BackgroundColor','w', 'String','(Filename=)',...
    'FontWeight','Bold', 'Enable','Inactive', 'ButtonDownFcn',@ClickFilename);

handles.ProgressText = uicontrol('Style','Text',...
    'Position',[20 100 360 20],'BackgroundColor','w', 'String','(Progress)');


% Input: Max Frame
handles.MaxFrameText = uicontrol('Style', 'Text', ...
    'Position', [20 80 150 20], 'BackgroundColor', 'w', ...
    'String', 'Max Frame:', 'HorizontalAlignment', 'left');
handles.MaxFrameEdit = uicontrol('Style', 'Edit', ...
    'Position', [180 80 200 20], 'BackgroundColor', 'w', ...
    'String', '8200'); % Default value for maxFrame

% Input: Add File Path
handles.AddFileText = uicontrol('Style', 'Text', ...
    'Position', [20 60 150 20], 'BackgroundColor', 'w', ...
    'String', 'Add File Path:', 'HorizontalAlignment', 'left');
handles.AddFileEdit = uicontrol('Style', 'Edit', ...
    'Position', [180 60 200 20], 'BackgroundColor', 'w', ...
    'String', ''); % Default value for AddFile

handles.RefImgLoaded = false;
handles.numPlanes = 1;  % default

% store handles in guidata
guidata(handles.fig, handles)



% callback function
function ClickStart(h, e)
    
    % retrieve guidata
    handles = guidata(h);
    % maxFrame=8200;
    % AddFile='E:\LuSLMOnlineTest\SL0838-Ai203\01142025\TSeries-01142025-1221-001';
    % Get maxFrame and AddFile from the GUI
    maxFrame = str2double(handles.MaxFrameEdit.String); % Convert to numeric
    AddFile = handles.AddFileEdit.String; % Get string input

    % Check if maxFrame is valid
    if isnan(maxFrame) || maxFrame <= 0
        error('Invalid value for Max Frame. Please enter a positive number.');
    end

    % Display inputs for debugging
    disp(['Max Frame: ', num2str(maxFrame)]);
    disp(['Add File: ', AddFile]);

% get ready for online registration
    DoRegistration = handles.DoRegistration.Value;  % get value, don't want to ever turn on or off mid acquisition.
    if DoRegistration
        if ~handles.RefImgLoaded
        LoadReferenceImage(h,e);
        end
    
         handles = guidata(h);
         numPlanes = handles.numPlanes;
         ops = handles.ops;
         refImg = handles.refImg;
    
    % initialise gpu. (is this needed?)
        for i = size(refImg,3)
%         gframe = single(gpuArray(refImg(:,:,1)));
           gframe = single(refImg(:,:,1));
        end
    % clear gframe?
    else
        numPlanes = 1;
    end

    % initialise PrairieLink
    pl = actxserver('PrairieLink.Application');
    pl.Connect();
    pl.SendScriptCommands('-DoNotWaitForScans');
    pl.SendScriptCommands('-LimitGSDMABufferSize true 100');
    pl.SendScriptCommands('-StreamRawData true 120');
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
    

    % write file header;  Commented by Lu Zhang, 2024-03-13;
    % fwrite(fileID, pixelsPerLine, 'uint16');
    % fwrite(fileID, linesPerFrame, 'uint16');

    % open binary file for writing
    if isempty(AddFile)
         if DoRegistration
            fileID = fopen([filePath '.bin'], 'wb');
            fileIDraw = fopen([filePath 'Raw.bin'], 'wb');
            shiftsAndCorrFileID = fopen([filePath '_ShiftsAndCorr.bin'],'wb');
        else
            fileID = fopen([filePath '.bin'], 'wb');
        end
    else
         if DoRegistration
            fileID = fopen([AddFile '.bin'], 'ab');
            fileIDraw = fopen([AddFile 'Raw.bin'], 'ab');
            shiftsAndCorrFileID = fopen([AddFile '_ShiftsAndCorr.bin'],'ab');
        else
            fileID = fopen([AddFile '.bin'], 'ab');
        end


    end



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
        subplot(1,2,1);
        Image1 = imagesc(zeros(linesPerFrame, pixelsPerLine));

        FrameCounter = title('');
        axis off; axis square; axis tight;
        subplot(1,2,2);
        Image2 = imagesc(zeros(linesPerFrame, pixelsPerLine));

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

%         disp([num2str(length(buffer)) ' of samples detected']);
%         Visulize samples recorded in each loop, only use for debugging!
%         Lu Zhang
%         figure(3);
%         subplot(3,1,1)
%         hold on;plot(icount,numSamplesRead,'r.')
%         subplot(3,1,2)
%         hold on;plot(icount,length(buffer),'r.')

%         subplot(3,1,3)
%         if exist('numWholeFramesGrabbed')
% %         hold on;plot(icount,numWholeFramesGrabbed,'b.')
%         icount=icount+1;
% 
%            hold on;plot(loopCounter,frameNum,'g.')
%            hold on;plot(loopCounter,framesCounter,'r.')
% 
% 
%         end
       
        % process the acquired frames (timer = ~5ms)
        if numWholeFramesGrabbed > 0
            for i = 1:numWholeFramesGrabbed
                if started == 0
                    started = 1;
                end
                 % get plane
                plane = mod(frameNum,numPlanes)+1;

                % get single frame
                frame = toProcess(((i-1)*totalSamplesPerFrame)+1:(i*totalSamplesPerFrame));

                % process the frame (C++ mex code)
                frame = PrairieLink_ProcessFrame(frame, samplesPerPixel, linesPerFrame, pixelsPerLine, flipEvenRows);

            % register frame HD 20180702
               if DoRegistration
%                 [regFrame,dv,cv] = return_offsets_phasecorr(single(gpuArray(frame)),ops{plane});
                  [regFrame,dv,cv] = return_offsets_phasecorr(single((frame)),ops{plane});
               
                % save processed frame and correlation values to file
                  fwrite(fileID, uint16(regFrame), 'uint16');
                  fwrite(fileIDraw, uint16(frame), 'uint16');

                  fwrite(shiftsAndCorrFileID, [gather(dv) gather(cv)], 'single');
               else
                 fwrite(fileID, frame, 'uint16');
               end


                % increment frame counter
                frameNum = frameNum + 1;


                % debugging: preview plot
                if preview
                    Image1.CData = [frame'];
                    Image2.CData = [ regFrame'];

                    FrameCounter.String = msg;
                    pause(0.00001);
                end

                if frameNum==maxFrame
                   break
                end
            end
        end

       % display progress
         if DoRegistration&exist('dv')
             msg = ['Frame: ' num2str(frameNum) ', Loop: ' num2str(loopCounter) '. Shifts: ' num2str(gather(dv(1)),'%.1f') ', ' num2str(gather(dv(2)),'%.1f') ];
         else
             msg = ['Frame: ' num2str(frameNum) ', Loop: ' num2str(loopCounter)];
         end
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

        if started && frameNum >= maxFrame
           running = 0;
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
    if DoRegistration
       fclose(shiftsAndCorrFileID);
       fclose(fileIDraw);

    end

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


function LoadReferenceImage(h, e)
% retrieve guidata
handles = guidata(h);

% select the image
[fileName,dirName] = uigetfile('*.tif','Select reference image/stack for registration','MultiSelect','on');
% cd(dirName)
if ~iscell(fileName)
    fileName = {fileName};
end

% build full paths
fullPath = [];
for i = 1:numel(fileName)
    fullPath{i} = [dirName filesep fileName{i}];
end

% set the gui label
fileNameList=[];
for i = 1:numel(fileName)
    fileNameList=[fileNameList fileName{i}];
end

handles.RefImgText.String = fileNameList;

% load image(s)
refImg = [];
ops = cell(numel(fullPath),1);
for i = 1:numel(fullPath)
    temp = imread(fullPath{i});
    refImg(:,:,i) = permute(temp,[2 1]);  % because PL data is different index order, but permuted for matlab when reading in.
    [ops{i,1}] = setup_registration_phasecorr(refImg(:,:,i),0);
end

% save to handles
handles.numPlanes = size(refImg,3);
handles.ops = ops;
handles.refImg = refImg;
handles.RefImgLoaded = true;

% set do registration to true (why else did you load the ref img?)
handles.DoRegistration.Value = true;

% store handles in guidata
guidata(handles.fig, handles)
end

