%%
% --- Lu Revision 09/6/2023, from original Ana's MotionCorrection2020MorePadding.m
% --- 09/14/2022 --- Just changed the upsamplign factor for testing
% --- Ana revision 06/15/2022 --- opts.DataRange = 'A2';
% --- Ana revision 05/30/2022 --- no change, just verification
% --- Ana revision 05/19/2022 --- comment line 352 (did not resolve issue)
% --- Ana revision 05/18/2022 --- removed excelSheetNum from input vars
% --- Ana revision 05/11/2022 --- xlsread replacement
% --- Ana revision 01/27/2022 --- our card reader chanel 2 pin broke, so
% our green chanel is now chanel 3, and therefore had to adapt the code.
% --- Ana revision 09/21/2021 --- (no change was made since final data
% processing, except for addition of more padding; all correct)

% This codes needs to be significantly improved
% At the moment is VERY inneficient for multiplane acquisition -> Suite2P

function [] = Meta2MatlabRW(tempCaDir, processedDir, excelfilename, rowsOfInterest)

% --- Get experiment folders under input directory ---
expFiles = dir(tempCaDir);
expFiles = expFiles(~strncmpi('.', {expFiles.name}, 1));

strncmpi('.', {expFiles.isdir}, 1);
% --- Read master excel file with experiment info ---
% [~,text,~] = xlsread(excelfilename,excelSheetNum);
opts = detectImportOptions(excelfilename); opts.DataRange = 'A2';
text = readmatrix(excelfilename, opts);

header = opts.VariableNames;
for i = 1:numel(header)
    if ~isempty(header{i})
        header{i} = strrep(header{i},' ','');
        param.(header{i}) = text(1:end,i);
    end
end
% header = text(1,:);
shift = 1;
% param = struct();
% for i = 1:numel(header)
%     if ~isempty(header{i})
%         header{i} = strrep(header{i},' ','');
%         param.(header{i}) = text(1+shift:end,i);
%     end
% end

clear opts header text

if isempty(rowsOfInterest)
    rowsOfInterest = 1:numel(param.Animal);
end

% --- Cycle through each experiment (i.e. animal, date, exp division) ---
% --- Treats each plane of each experiment independently ---
for m = 1:numel(rowsOfInterest)
    
    k = rowsOfInterest(m)-shift;
    
    tempAnimal=param.Animal{k};
    i1=strfind(tempAnimal,'-');
    if ~isempty(i1)
       tempAnimal=tempAnimal(1:i1-1);
    end
    
    currExp = expFiles(contains({expFiles.name}, tempAnimal) == 1 & contains({expFiles.name}, param.Date{k}) == 1);

   for ii = 1:length(currExp)
       if currExp(ii).isdir
          currExp = currExp(ii);
          break
       end
    end


    pathGet = [currExp.folder '\' currExp.name];
    
    strg = strfind(pathGet,'\');

    %% Define the saving file name
    % if isfield(param,'ExpDivision') == 1
    %     experiment = [pathGet(strg(end)+1:end) '_' param.ExpDivision{k,1}];
    % else
    %     experiment = pathGet(strg(end)+1:end);
    % end
    if isfield(param,'ExpDivision') == 1
        experiment = [param.Date{k,1} '_' param.Animal{k,1} '_' param.Date_2{k,1} '_' param.ExpDivision{k,1}];
    else
        experiment = [param.Date{k,1} '_' param.Animal{k,1} '_' param.Date_2{k,1}];
    end

    %% 
    clear strg
    
    disp(['---Preprocessing ' experiment '---'])
    tic
    allFiles = dir([pathGet '\']);
    allFiles = allFiles(~strncmpi('.', {allFiles.name}, 1));
    caTrialsTemp = allFiles(strncmp('TSeries', {allFiles.name}, numel('TSeries')));
    caTrialsInd = str2num(param.Tseries{k,1});
    clear caTrials
    for i = 1:numel(caTrialsInd)
        caTrials(i) = caTrialsTemp(endsWith({caTrialsTemp.name}, sprintfc('%0.3d',caTrialsInd(i))));
    end
    
    clear allFiles CaTrialsTemp
    
    % caTrialsSpont = caTrials(endsWith({caTrials.name}, sprintfc('%0.3d',str2num(param.TseriesSpont{k,1}))));
    % caTrialsEvoked = caTrials(endsWith({caTrials.name}, sprintfc('%0.3d',str2num(param.TseriesEvoked{k,1}))));
    
    %% Check if there is Sound trial, I assume it is the CS-reward?
    caTrialsSound = caTrials(endsWith({caTrials.name}, sprintfc('%0.3d',str2num(param.TseriesSound{k,1}))));
    

    %% Check if there is open opto light control?
    if isfield(param, 'TseriesOpto')
        caTrialsLED = caTrials(endsWith({caTrials.name}, sprintfc('%0.3d',str2num(param.TseriesOpto{k,1}))));
    else
        caTrialsLED = [];
    end
    

    %% Check if there is close-loop light control
    caTrialsLEDCL = caTrials(endsWith({caTrials.name}, sprintfc('%0.3d',str2num(param.TseriesOptoCL{k,1}))));
    
    %% Check if there is red fluorescence data
    caTrialsRed = caTrials(endsWith({caTrials.name}, sprintfc('%0.3d',str2num(param.TrialsWithRed{k,1}))));
    caTrialsRedMask = zeros(length(caTrials),1);
    for i = 1:length(caTrialsRed)
        index = strcmp({caTrials.name}, caTrialsRed(i).name) == 1;
        caTrialsRedMask(index) = 1;
    end
    
    clear caTrialsRed index
    
%% Read the Mega data
    % --- Get CaRec, VoltRec and VoltOutput metadata ---
    % --- I keep it in memory and save it later for each plane (bad) ---
    vRec = [];
    vOut = [];
    infoAllTrials = [];
    for i = 1:numel(caTrialsInd)
        caTrialCurr = caTrials(endsWith({caTrials.name}, sprintfc('%0.3d',caTrialsInd(i))));
        fieldname = char(sprintfc('Trial%0.3d',caTrialsInd(i)));
        
        % --- Read calcium recordings metadata ---
        infoAllTrials.(fieldname) = xml2struct([pathGet '\' caTrialCurr.name '\' caTrialCurr.name '.xml']);
        % --- Read voltage recordings metadata ---
        vFile = dir([pathGet '\' caTrialCurr.name '\*VoltageRecording*.csv']);
        vRec.(fieldname) = csvread([pathGet '\' caTrialCurr.name '\' vFile.name],1,0);
        
        % --- Read voltage output metadata ---
        % --- Use this data to confirm master excel table ---
        vOutFile = dir([pathGet '\' caTrialCurr.name '\*VoltageOutput*.xml']);
        vOutTemp = xml2struct([vOutFile.folder '\' vOutFile.name]);
        caTrials(i).LED = vOutTemp.Experiment.Waveform{1,2}.Enabled.Text;
        caTrials(i).Evoked = vOutTemp.Experiment.Waveform{1,3}.Enabled.Text;
        
        caTrials(i).LED_Delay = vOutTemp.Experiment.Waveform{1,2}.WaveformComponent_PulseTrain.FirstPulseDelay.Text;
        caTrials(i).LED_Num = vOutTemp.Experiment.Waveform{1,2}.WaveformComponent_PulseTrain.PulseCount.Text;
        caTrials(i).LED_Width = vOutTemp.Experiment.Waveform{1,2}.WaveformComponent_PulseTrain.PulseWidth.Text;
        caTrials(i).LED_ISI = vOutTemp.Experiment.Waveform{1,2}.WaveformComponent_PulseTrain.PulseSpacing.Text;
        
        caTrials(i).Evoked_Delay = vOutTemp.Experiment.Waveform{1,3}.WaveformComponent_PulseTrain.FirstPulseDelay.Text;
        caTrials(i).Evoked_Num = vOutTemp.Experiment.Waveform{1,3}.WaveformComponent_PulseTrain.PulseCount.Text;
        caTrials(i).Evoked_Width = vOutTemp.Experiment.Waveform{1,3}.WaveformComponent_PulseTrain.PulseWidth.Text;
        caTrials(i).Evoked_ISI = vOutTemp.Experiment.Waveform{1,3}.WaveformComponent_PulseTrain.PulseSpacing.Text;
        vOut.(fieldname) = vOutTemp;
    end
    
    clear vFile vOutFile vOutTemp
    
    fieldname = char(sprintfc('Trial%0.3d',str2num(param.TrialForReg{k,1})));
    
    % We assume that numPlanes, num_x and num_y are constant across the
    % experiment (these parameters MUST be constant across trials)
    numPlanes = 1;
    if size(infoAllTrials.(fieldname).PVScan.Sequence,2) > 1
        numPlanes = size(infoAllTrials.(fieldname).PVScan.Sequence{1, 1}.Frame,2);
    end
    num_x = str2num(infoAllTrials.(fieldname).PVScan.PVStateShard.PVStateValue{1, 10}.Attributes.value);  % Lines per frame
    num_y = str2num(infoAllTrials.(fieldname).PVScan.PVStateShard.PVStateValue{1, 18}.Attributes.value);  % Pixels per line
    
    % --- Load frames for generation of reference images ChX ---
    % --- We typically use one trial, either Ch1 (red) or Ch2/Ch3 (green) ---
    disp('---Load caTrialReg---')
    caTrialReg = caTrials(endsWith({caTrials.name}, sprintfc('%0.3d',str2num(param.TrialForReg{k,1}))));
    [framesTrialReg] = readTifsAna2019([pathGet '\' caTrialReg.name], num_x, num_y, ['Ch' param.ChForReg{k,1}]);
    
    % --- Median filter data to remove shot noise ---
    % framesAllFilt = [];
    % for i = 1:num_z
    %     [allFramesFilt(:,:,i)] = medfilt2(framesAll(:,:,i));
    % end
    
    % --- Treat each plane as an independent experiment ---
    % --- Then, for each plane, cycle through trials ---
    for selectedPlane = 1:numPlanes % str2num(param.PlaneN{k,1})
        
        % [pathSave] = [processedDir experiment '_Plane' sprintf('%02i',selectedPlane) '_Analysed'];
        [pathSave] = [processedDir];
        
        if ~exist(pathSave,'dir')
            mkdir(pathSave)
        end
        
        if exist(pathSave,'dir')
            filesToDel = dir([pathSave '\*.tif']);
            for d = 1:length(filesToDel)
                delete([pathSave '\' filesToDel(d).name]);
            end
        end
        clear d
        
        % --- Get frame start time stamp and frame period ---
        % numPlanes = []; num_x = []; num_y = [];
        timeStampCa = [];
        lastFrame = [];
        lastFrame(1) = 0;
        for i = 1:numel(caTrialsInd)
            
            caTrialCurr = caTrials(endsWith({caTrials.name}, sprintfc('%0.3d',caTrialsInd(i))));
            fieldname = char(sprintfc('Trial%0.3d',caTrialsInd(i)));
            
            info = infoAllTrials.(fieldname);
            
            % numPlanes(i) = size(info.PVScan.Sequence,2);
            % num_x(i) = str2num(info.PVScan.PVStateShard.PVStateValue{1, 10}.Attributes.value);  % Lines per frame
            % num_y(i) = str2num(info.PVScan.PVStateShard.PVStateValue{1, 18}.Attributes.value);  % Pixels per line
            
            timeStampCaTemp = [];
            if numPlanes == 1
                testError = size(info.PVScan.Sequence,2);
                if testError == 1
                    for ii = 1:numel(info.PVScan.Sequence.Frame) % Number of frames
                        timeStampCaTemp(ii,1) = str2num(info.PVScan.Sequence.Frame{1, ii}.Attributes.relativeTime);
                        timeStampCaTemp(ii,2) = str2num(info.PVScan.Sequence.Frame{1, ii}.Attributes.absoluteTime); % in sec
                        lastFrame(i+1,:) = length(timeStampCa);
                    end
                end
                if testError >1
                    for ii = 1:numel(info.PVScan.Sequence{1,1}.Frame) % Number of frames
                        timeStampCaTemp(ii,1) = str2num(info.PVScan.Sequence{1,1}.Frame{1, ii}.Attributes.relativeTime);
                        timeStampCaTemp(ii,2) = str2num(info.PVScan.Sequence{1,1}.Frame{1, ii}.Attributes.absoluteTime); % in sec
                        lastFrame(i+1,:) = length(timeStampCa);
                    end
                end
            else
                for ii = 1:numel(info.PVScan.Sequence) % Number of cycles, frame "selectedPlane" within each cycle
                    timeStampCaTemp(ii,1) = str2num(info.PVScan.Sequence{1, ii}.Frame{1, selectedPlane}.Attributes.relativeTime);
                    timeStampCaTemp(ii,2) = str2num(info.PVScan.Sequence{1, ii}.Frame{1, selectedPlane}.Attributes.absoluteTime); % in sec
                    if size(info.PVScan.Sequence{1, ii}.Frame{1, selectedPlane}.PVStateShard.PVStateValue,2)>1
                        timeStampCaTemp(ii,3) = str2num(info.PVScan.Sequence{1, ii}.Frame{1, selectedPlane}.PVStateShard.PVStateValue{1, 1}.Attributes.value);
                    else
                        timeStampCaTemp(ii,3) = str2num(info.PVScan.Sequence{1, ii}.Frame{1, selectedPlane}.PVStateShard.PVStateValue.Attributes.value); % I found this issue in the metadata in that the fields are not consistent (in some experiments) -> why?
                    end
                end
            end
            timeStampCa = [timeStampCa; timeStampCaTemp];
            lastFrame(i+1,:) = length(timeStampCa);
            num_z = size(timeStampCa,1);
        end
        
        % --- Save variables I ---
        save([pathSave '\' experiment '_Plane' sprintf('%02i',selectedPlane) '_vars'], 'experiment', 'pathGet', 'pathSave', 'caTrials', 'caTrialsSound', 'caTrialsLED', 'caTrialsLEDCL', ...
            'num_x','num_y','num_z','lastFrame','timeStampCa','vRec','vOut') %,'infoAllTrials')

        toc
    end
end