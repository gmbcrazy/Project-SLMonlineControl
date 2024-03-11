%%
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

function [] = MotionCorrection2020MorePadding_updateHari(tempCaDir, processedDir, excelfilename, rowsOfInterest)

% --- Get experiment folders under input directory ---
expFiles = dir(tempCaDir);
expFiles = expFiles(~strncmpi('.', {expFiles.name}, 1));

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
clear i
if isempty(rowsOfInterest)
    rowsOfInterest = 1:numel(param.Animal);
end

% --- Cycle through each experiment (i.e. animal, date, exp division) ---
% --- Treats each plane of each experiment independently ---
for m = 1:numel(rowsOfInterest)
    
    k = rowsOfInterest(m)-shift;
    
    currExp = expFiles(contains({expFiles.name}, param.Animal{k}) == 1 & contains({expFiles.name}, param.Date{k}) == 1);
    pathGet = [currExp.folder '\' currExp.name];
    
    strg = strfind(pathGet,'\');
    if isfield(param,'ExpDivision') == 1
        experiment = [pathGet(strg(end)+1:end) '_' param.ExpDivision{k,1}];
    else
        experiment = pathGet(strg(end)+1:end);
    end
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
%%%%comments by Hari 82-100    
% % %     caTrialsSound = caTrials(endsWith({caTrials.name}, sprintfc('%0.3d',str2num(param.TseriesSound{k,1}))));
% % %     
% % %     if isfield(param, 'TseriesOpto')
% % %         caTrialsLED = caTrials(endsWith({caTrials.name}, sprintfc('%0.3d',str2num(param.TseriesOpto{k,1}))));
% % %     else
% % %         caTrialsLED = [];
% % %     end
% % %     
% % %     caTrialsLEDCL = caTrials(endsWith({caTrials.name}, sprintfc('%0.3d',str2num(param.TseriesOptoCL{k,1}))));
% % %     
% % %     caTrialsRed = caTrials(endsWith({caTrials.name}, sprintfc('%0.3d',str2num(param.TrialsWithRed{k,1}))));
% % %     caTrialsRedMask = zeros(length(caTrials),1);
% % %     for i = 1:length(caTrialsRed)
% % %         index = strcmp({caTrials.name}, caTrialsRed(i).name) == 1;
% % %         caTrialsRedMask(index) = 1;
% % %     end
% % %     
% % %     clear caTrialsRed index
    
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
        
       % caTrials(i).LED_Delay = vOutTemp.Experiment.Waveform{1,2}.WaveformComponent_PulseTrain.FirstPulseDelay.Text;
        %caTrials(i).LED_Num = vOutTemp.Experiment.Waveform{1,2}.WaveformComponent_PulseTrain.PulseCount.Text;
        %caTrials(i).LED_Width = vOutTemp.Experiment.Waveform{1,2}.WaveformComponent_PulseTrain.PulseWidth.Text;
        %caTrials(i).LED_ISI = vOutTemp.Experiment.Waveform{1,2}.WaveformComponent_PulseTrain.PulseSpacing.Text;
        
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
        
        [pathSave] = [processedDir experiment '_Plane' sprintf('%02i',selectedPlane) '_Analysed'];
        
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
        %%change by Hari 232 to 233
        save([pathSave '\' experiment '_Plane' sprintf('%02i',selectedPlane) '_vars'], 'experiment', 'pathGet', 'pathSave', 'caTrials', ...
            'num_x','num_y','num_z','lastFrame','timeStampCa','vRec','vOut') %,'infoAllTrials')
        %save([pathSave '\' experiment '_Plane' sprintf('%02i',selectedPlane) '_vars'], 'experiment', 'pathGet', 'pathSave', 'caTrials', 'caTrialsSound', 'caTrialsLED', 'caTrialsLEDCL', ...
            %'num_x','num_y','num_z','lastFrame','timeStampCa','vRec','vOut') %,'infoAllTrials')
        
        % --- Create reference image per plane ChX ---
        framesPlane = framesTrialReg(:,:,selectedPlane:numPlanes:size(framesTrialReg,3)); % framesTrialReg(:,:,selectedPlane:str2num(param.PlaneN{k,1}):size(framesTrialReg,3));
        avFrame = mean(framesPlane,3);
        avFrameFT = fft2(avFrame);
        imwrite (uint16(avFrame), [pathSave '\' experiment '_Plane' sprintf('%02i',selectedPlane) '_avFrame.tif'])
        
        clear framesPlane avFrame
        
        % --- Registration parameters ---
        usfac = 30; % 20; % upsampling factor (integer) for dftregistration
        pd = 80; % 40; % padding for shifting frames (motion correction)
        hpd = pd/2; % half padding
        
        offsetsCa = []; % should allocate
        
        % --- Cycle through trials for each plane ---
        % --- Registration and alignment of Ch1 alone or Ch1 and Ch2/Ch3 ---
        disp(['---Motion Correction_Plane' num2str(selectedPlane) '---'])
        for i = 1:numel(caTrialsInd)
            
            caTrialCurr = caTrials(endsWith({caTrials.name}, sprintfc('%0.3d',caTrialsInd(i))));
            fieldname = char(sprintfc('Trial%0.3d',caTrialsInd(i)));
            
            % --- Load trial ChX ---
            [framesTrial] = readTifsAna2019([pathGet '\' caTrialCurr.name], num_x, num_y, ['Ch' param.ChForReg{k,1}]);
            framesTrial = framesTrial(:,:,selectedPlane:str2num(param.PlaneN{k,1}):size(framesTrial,3));
            
            % --- Register ChX ---
            for ii = 1:size(framesTrial,3)
                % The following part was for optogenetics, but in general is
                % removing frames in noisy records -> Do it a different way
                % if mean(mean(framesTrial(:,:,ii)))<200
                %   offsetsCa.(fieldname)(ii,:) = [0 0 0 0];
                % else
                % offsetsCa.(fieldname)(ii,:) = zeros(4500,4);
                offsetsCa.(fieldname)(ii,:) = dftregistration(avFrameFT, fft2(framesTrial(:,:,ii)), usfac);  % output of dftregistration [error,diffphase,net_row_shift,net_col_shift]
            end
            
            % --- Align ChX---
            temp = uint16(zeros(num_x+pd,num_y+pd,size(framesTrial,3)));
            for ii = 1:size(framesTrial,3)
                sx = round(offsetsCa.(fieldname)(ii,3));
                sy = round(offsetsCa.(fieldname)(ii,4));
                % Check if shifts are within the allowed range (1:num_x+pd/1:num_y+pd)
                if hpd+sx > 0 && hpd+sx+num_x-1 <= num_x+pd
                    if hpd+sy > 0 && hpd+sy+num_y-1 <= num_y+pd
                        temp(sx+hpd:sx+hpd+num_x-1, sy+hpd:sy+hpd+num_y-1, ii)...
                            = framesTrial(:,:,ii);
                    else
                        temp (:,:,ii) = 0;
                    end
                else
                    temp (:,:,ii) = 0;
                end
            end
            
            % --- Clear uncorrected matrices ---
            clear framesTrial
            
            % --- Remove padding/edges ---
            temp = temp(hpd:num_x+hpd-1, hpd:num_y+hpd-1, :);
            
            % --- Create aligned trial ChX average and std ---
            avTemp = mean(temp,3);
            imwrite (uint16(avTemp), [pathSave '\' experiment '_Plane' sprintf('%02i',selectedPlane) '_Ch' param.ChForReg{k,1} '_AV.tif'], 'WriteMode', 'append')
            stdTemp = std(double(temp), 0, 3);
            imwrite (uint16(stdTemp), [pathSave '\' experiment '_Plane' sprintf('%02i',selectedPlane) '_Ch' param.ChForReg{k,1} '_STD.tif'], 'WriteMode', 'append')
            
            if str2num(param.ChForReg{k,1}) == 2 || str2num(param.ChForReg{k,1}) == 3 % If ChX is Ch2/Ch3 (green)
                
                % --- Save aligned images Ch2/Ch3 ---
                % writeTifChunks2019(temp, [pathSave '\' experiment '_Plane' sprintf('%02i',selectedPlane) '_' fieldname '_MC'], num_z/(length(lastFrame)-1)/5);
                
                % --- Clear uncorrected matrices II ---
                clear temp avTemp stdTemp
              %% commented by Hari from 313-350 
                % --- Align Ch1 in trials that have it ---
%                 caTrialsRedMask=[];
%                 if caTrialsRedMask(i) == 1
%                     
%                     % --- Load and align Ch1 (red) ---
%                     framesTrial1 = readTifsAna2019([pathGet '\' caTrialCurr.name], num_x, num_y, 'Ch1');
%                     framesTrial1 = framesTrial1(:,:,selectedPlane:str2num(param.PlaneN{k,1}):size(framesTrial1,3));
%                     
%                     temp1 = uint16(zeros(num_x+pd,num_y+pd,size(framesTrial1,3)));
%                     
%                     for ii = 1:size(framesTrial1,3)
%                         sx = round(offsetsCa.(fieldname)(ii,3));
%                         sy = round(offsetsCa.(fieldname)(ii,4));
%                         % Check if shifts are within the allowed range (1:num_x+pd/1:num_y+pd)
%                         if hpd+sx > 0 && hpd+sx+num_x-1 <= num_x+pd
%                             if hpd+sy > 0 && hpd+sy+num_y-1 <= num_y+pd
%                                 temp(sx+hpd:sx+hpd+num_x-1, sy+hpd:sy+hpd+num_y-1, ii)...
%                                     = framesTrial1(:,:,ii);
%                             else
%                                 temp1 (:,:,ii) = 0;
%                             end
%                         else
%                             temp1 (:,:,ii) = 0;
%                         end
%                     end
%                     % --- Clear uncorrected matrices I ---
%                     clear framesTrial1
%                     
%                     % --- Remove padding/edges ---
%                     temp1 = temp1(hpd:num_x+hpd-1, hpd:num_y+hpd-1, :);
%                     
%                     % --- Create trial average and std ---
%                     avTemp1 = mean(temp1,3);
%                     imwrite (uint16(avTemp1), [pathSave '\' experiment '_Plane' sprintf('%02i',selectedPlane) '_Ch1_AV.tif'], 'WriteMode', 'append')
%                     
%                     % --- Clear uncorrected matrices II ---
%                     clear temp1 avTemp1
%                 end
            else % If ChX is Ch1 (red)
                
                % --- Clear uncorrected matrices II ---
                clear temp avTemp stdTemp
                
                % --- Load and align Ch2/Ch3 (green) ---
                framesTrial1 = readTifsAna2019([pathGet '\' caTrialCurr.name], num_x, num_y, ['Ch' param.ChGreen{k,1}]); % could find it automatically
                framesTrial1 = framesTrial1(:,:,selectedPlane:str2num(param.PlaneN{k,1}):size(framesTrial1,3));
                
                temp1 = uint16(zeros(num_x+pd,num_y+pd,size(framesTrial1,3)));
                
                for ii = 1:size(framesTrial1,3)
                    sx = round(offsetsCa.(fieldname)(ii,3));
                    sy = round(offsetsCa.(fieldname)(ii,4));
                    % Check if shifts are within the allowed range (1:num_x+pd/1:num_y+pd)
                    if hpd+sx > 0 && hpd+sx+num_x-1 <= num_x+pd
                        if hpd+sy > 0 && hpd+sy+num_y-1 <= num_y+pd
                            temp1(sx+hpd:sx+hpd+num_x-1, sy+hpd:sy+hpd+num_y-1, ii)...
                                = framesTrial1(:,:,ii);
                        else
                            temp1 (:,:,ii) = 0;
                        end
                    else
                        temp1 (:,:,ii) = 0;
                    end
                end
                
                % --- Clear uncorrected matrices I ---
                clear framesTrial1
                
                % --- Remove padding/edges ---
                temp1 = temp1(hpd:num_x+hpd-1, hpd:num_y+hpd-1, :);
                
                % --- Create trial average and std ---
                avTemp1 = mean(temp1,3);
                imwrite (uint16(avTemp1), [pathSave '\' experiment '_Plane' sprintf('%02i',selectedPlane) '_Ch2_AV.tif'], 'WriteMode', 'append')
                stdTemp1 = std(double(temp1), 0, 3);
                imwrite (uint16(stdTemp1), [pathSave '\' experiment '_Plane' sprintf('%02i',selectedPlane) '_Ch2_STD.tif'], 'WriteMode', 'append')
                
                % --- Save aligned images (just for the moment) ---
                % writeTifChunks2019(temp1, [pathSave '\' experiment '_Plane' sprintf('%02i',selectedPlane) '_' fieldname '_MC'], num_z/(length(lastFrame)-1)/5);
                
                % --- Clear uncorrected matrices II ---
                clear temp1 avTemp1 stdTemp1
            end
        end
        
        % --- Save variables II ---
        save([pathSave '\' experiment '_Plane' sprintf('%02i',selectedPlane) '_vars'], 'offsetsCa', '-append');
        toc
    end
end