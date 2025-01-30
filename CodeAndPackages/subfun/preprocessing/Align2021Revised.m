%%

% --- Ana last revision 9/26/2021 --- (I worked on it 9/29 but made no
% change)

% Takes in the Ca and Whisk .mat files, aligns, crops, and resamples raw
% variables to generate working variables for analysis

% Prelast varin = 0 or 1, vRec contains LFP or not, respectively
% Should create an option to read it directly from the excel master
% In future improve the structure

% Eventually, remove legacy part kept for motion correction without added
% column 3 in timeStampCa (which is frame period)

function [] = Align2021Revised(currExpVars, deleteColumn2)

% --- Load vars ---
load([currExpVars.folder '\' currExpVars.name])
%load([currExpVarsWhisk.folder '\' currExpVarsWhisk.name])%%by Hari

% --- Gather information about each trial identity ---

% caTrials -> All trials taken into analysis
% caTrialSound -> Manually saved information needed to categorize trials
% All remaining information is automatically read from vOut
% I.e., we can read from the metadata if the whisker stimulation is on
% (stim), but we can only determine by the actual protocol if the whisker
% stimulator was (stim) or not (sound) touching the whiskers

caTrialsMask = [];
caTrialsMask(1:length(caTrials),1) = 0;
% 0 = spontTrials; 1 = stimTrials; 2 = soundTrials
for t = 1:length(caTrials)
    if strcmp((caTrials(t).Evoked),'true') == 1
        caTrialsMask(t,1) = 1;
    end
end
%%commented by Hari 40-45
% if ~isempty(caTrialsSound)
%     for t = 1:length(caTrialsSound)
%         index = find(strcmp({caTrials.name}, caTrialsSound(t).name) == 1);
%         caTrialsMask(index,1) = 2;
%     end
% end

clear index

caTrialsMaskLED = [];
caTrialsMaskLED(1:length(caTrials),1) = 0;
if exist('caTrialsLED','var')==1 && ~isempty(caTrialsLED) % I added if exists because in most runs, I read it automatically from metadata, but had to create this field for the neuromod package
    for t = 1:length(caTrialsLED)
        index = find(strcmp({caTrials.name}, caTrialsLED(t).name) == 1);
        caTrialsMaskLED(index,1) = 1;
    end
else
    for t = 1:length(caTrials)
        if strcmp((caTrials(t).LED),'true') == 1
            caTrialsMaskLED(t,1) = 1;
        end
    end
end
%%commented by Hari 64-72
% caTrialsMaskLEDCL = [];
% caTrialsMaskLEDCL(1:length(caTrials),1) = 0;
% if ~isempty(caTrialsLEDCL)
%     for t = 1:length(caTrialsLEDCL)
%         index = find(strcmp({caTrials.name}, caTrialsLEDCL(t).name) == 1);
%         caTrialsMaskLEDCL(index,1) = 1;
%     end
% end
% clear index

% ---
lastFrame(:,2:3) = 0;

% ---
%commented by Hari
% if size(posi,1)>1
%     whiskTrial1 = whiskTrial(:,:,1);
%     whiskTrial2 = whiskTrial(:,:,2);
% else
%     whiskTrial1 = whiskTrial;
% end
% clear whiskTrial

% ---
spks = spks;
fDeltaFoF = [];
fSpks = [];
fWhisk1Orig = [];
fWhisk1 = [];
fWhisk2 = [];
fSpeed = [];
fStim = [];
fSound = [];
fLED = [];
fLEDCL = [];
fTrial = [];

% --- Process trial by trial (t) ---
field = fieldnames(vRec);

for t = 1:length(field)
    
    fieldname = char(field(t));
    vRecTemp = vRec.(fieldname);
    
    % Delete column 2
    if deleteColumn2 == 1
        vRecTemp = vRecTemp(:,[1 3 4]);
    end
    
    % Add more working columns
    if size(vRecTemp,2) == 3 % Time, Speed, Cam
        vRecTemp(:,4:12) = NaN; vRecTemp(:,[7 10 11]) = 0;
    end
    if size(vRecTemp,2) == 4 % Time, Speed, Cam, LED
        vRecTemp(:,5:12) = NaN; vRecTemp(:,[7 10 11]) = 0;
    end
    
    % Clear column 4 if LEDCL was not on
    %(sometimes I forgot to disconnect the feedback cable)
    %%commnted by Hari 125-127
%     if caTrialsMaskLEDCL(t) == 0
%         vRecTemp(:,4) = 0;
%     end
    
    sampFreq = 1/diff(vRecTemp(1:2,1))*1000;
    
    % Cam rising edges matched to nFramesWhisk per trial
    %%commented by Hari133-149
%     fallingEdges = LocalMinima(diff(vRecTemp(:,3)),1,-1)+1;
%     fW = numel(find(whiskTrial1(:,t)>0))+1; % Total number of whisk frames actually recorded per trial; +1 because frame 1 is NaN
%     if numel(fallingEdges) > fW % It is possible because the cam stays live after recording is stopped
%         vRecTemp(fallingEdges(1:fW),5) = 1:fW;
%         vRecTemp(fallingEdges(1:fW),6) = whiskTrial1(1:fW,t);
%         vRecTemp(fallingEdges(1:fW),8) = whiskTrial1(1:fW,t);
%         if exist('whiskTrial2', 'var') == 1
%             vRecTemp(fallingEdges(1:fW),9) = whiskTrial2(1:fW,t); % First value of whisk for each trial is an artifact
%         end
%     else
%         vRecTemp(fallingEdges,5) = 1:numel(fallingEdges);
%         vRecTemp(fallingEdges,6) = whiskTrial1(1:numel(fallingEdges),t);
%         vRecTemp(fallingEdges,8) = whiskTrial1(1:numel(fallingEdges),t);
%         if exist('whiskTrial2', 'var') == 1
%             vRecTemp(fallingEdges,9) = whiskTrial2(1:numel(fallingEdges),t);
%         end
%     end
    
    % Check if the trial has stim or sound
    %%C By Hari 153-199
%     if caTrialsMask(t) > 0
%         pulseDelay   = str2num(vOut.(fieldname).Experiment.Waveform{1,3}.WaveformComponent_PulseTrain.FirstPulseDelay.Text)*(sampFreq/1000); % ms
%         pulseCount   = str2num(vOut.(fieldname).Experiment.Waveform{1,3}.WaveformComponent_PulseTrain.PulseCount.Text);
%         pulseWidth   = str2num(vOut.(fieldname).Experiment.Waveform{1,3}.WaveformComponent_PulseTrain.PulseWidth.Text)*(sampFreq/1000); % ms
%         pulseSpacing = str2num(vOut.(fieldname).Experiment.Waveform{1,3}.WaveformComponent_PulseTrain.PulseSpacing.Text)*(sampFreq/1000); % ms
%         
%         stimTrace = [];
%         stimTrace(1:size(vRecTemp,1),1) = 0;
%         stimTrace(pulseDelay+1 : pulseWidth+pulseSpacing : size(vRecTemp,1), 1) = 1;
%         
%         positions = find(stimTrace>0);
%         if numel(positions) > pulseCount
%             stimTrace(positions(pulseCount+1):end) = 0;
%             positions = find(stimTrace>0);
%         end 
%         for p = 1:numel(positions)
%             stimTrace(positions(p) : positions(p)+pulseWidth,1) = 1;
%         end
%         stimTrace = stimTrace(1:size(vRecTemp,1));
%         
%         % If the trial has stim, remove the artifact from the whiskTrial1
%         if caTrialsMask(t) == 1
%             vRecTemp(:,7) = stimTrace;
%             
%             for p = 1:numel(positions)
%                 stimTrace(positions(p) : positions(p)+2.5*pulseWidth,1) = 1;
%             end
%             stimTrace = stimTrace(1:size(vRecTemp,1));
% 
%             x = vRecTemp(:,5); x1 = x(x>0);
%             x2 = x; x2(stimTrace==1) = NaN; x3 = x2(x2>0);
%             y = whiskTrial1(x3,t);
%             yClean = interp1(x3, y, x1, 'linear');
%             
%             vRecTemp(find(vRecTemp(:,5)>0),8) = yClean(x1);
%             
%             if exist('whiskTrial2', 'var') == 1
%                 y2 = whiskTrial2(x3,t);
%                 yClean2 = interp1(x3, y2, x1, 'linear');
%                 vRecTemp(find(vRecTemp(:,5)>0),9) = yClean2(x1);
%             end
%         end
%         
%         % If the trial has sound, keep the original WhiskTrial1
%         if caTrialsMask(t) == 2
%             vRecTemp(:,10) = stimTrace;
%         end
%     end
   
    % Check if the trial has LED open-loop
    %%C By Hari 204-224
%     if strcmp((caTrials(t).LED),'true') == 1
%         pulseDelayLED   = str2num(vOut.(fieldname).Experiment.Waveform{1,2}.WaveformComponent_PulseTrain.FirstPulseDelay.Text)*(sampFreq/1000); % ms
%         pulseCountLED   = str2num(vOut.(fieldname).Experiment.Waveform{1,2}.WaveformComponent_PulseTrain.PulseCount.Text);
%         pulseWidthLED   = str2num(vOut.(fieldname).Experiment.Waveform{1,2}.WaveformComponent_PulseTrain.PulseWidth.Text)*(sampFreq/1000); % ms
%         pulseSpacingLED = str2num(vOut.(fieldname).Experiment.Waveform{1,2}.WaveformComponent_PulseTrain.PulseSpacing.Text)*(sampFreq/1000); % ms
%         
%         LEDTrace = [];
%         LEDTrace(1:size(vRecTemp,1),1) = 0;
%         LEDTrace(pulseDelayLED+1 : pulseWidthLED+pulseSpacingLED : size(vRecTemp,1), 1) = 1;
%         
%         positionsLED = find(LEDTrace>0);
%         if numel(positionsLED) > pulseCountLED
%             LEDTrace(positionsLED(pulseCountLED+1):end) = 0;
%             positionsLED = find(LEDTrace>0);
%         end
%         for p = 1:numel(positionsLED)
%             LEDTrace(positionsLED(p) : positionsLED(p)+pulseWidthLED,1) = 1;
%         end
%         LEDTrace = LEDTrace(1:size(vRecTemp,1));
%         
%         vRecTemp(:,11) = LEDTrace;
%     end
        
    % Add calcium time stamps and frame number
    timeStampCaTrial = timeStampCa(lastFrame(t)+1:lastFrame(t+1),1)*1000;
    timeStampCaTrial(length(timeStampCaTrial)+1) = timeStampCa(lastFrame(t+1))*1000+mean(diff(timeStampCaTrial(:,1)));
    
    addTime = timeStampCaTrial(2:end);
    addTime(:,2:12) = NaN;
    vRecTemp = sortrows([addTime; vRecTemp]);
    clear addTime
    
    fCa = lastFrame(t+1)-(lastFrame(t)+1)+1;
    frame = lastFrame(t)+1;
    [~,idx] = intersect(vRecTemp(:,1),timeStampCaTrial,'stable');
    % if fR > 25 NEED TO WORK ON THE MULTIPLANE PART
        for n = 1:fCa
            vRecTemp(idx(n):idx(n+1)-1,12) = frame;
            frame = frame + 1;
        end
    % else % NEED TO WORK ON THE MULTIPLANE PART
    %    for n = 1:length(timeStampCaTrial)
    %        vRecNew.(fieldname)(timeStampCaTrial(n,1)+1:(timeStampCaTrial(n,1)+1+timeStampCaTrial(n,3)),6) = frame;
    %        frame = frame + 1;
    %    end
    % end
    
    % Establish usable trial time
    startCa = min(find(vRecTemp(:,12) == lastFrame(t)+1));
    endCa = max(find(vRecTemp(:,12) == lastFrame(t+1)));
    %startWhisk = find(vRecTemp(:,5) == 1);C By Hari
    %endWhisk = find(vRecTemp(:,5) == max(vRecTemp(:,5)));C By Hari
    %C By Hari 257-261
%     if startWhisk > startCa % means that Whisk aquisition starts after Ca aquisition -> thus is limiting factor; typical stituation except in volume imaging
%         newStart = startWhisk;
%     else
%         newStart = startCa;
%     end
    %By Hari 263-271
%     if endWhisk > endCa % means that Ca aquisition ended before Whisk aquisition -> thus is limiting factor; typical stituation except in some first experiments in which I recorded only 37500 Whisk frames
%         % if fR > 25 NEED TO WORK ON THE MULTIPLANE PART
%             newEnd = endCa;
%         % else
%         %    newEnd = endCa + timeStampCaTrial(end,3); NEED TO WORK ON THE MULTIPLANE PART
%         % end
%     else
%         newEnd = endWhisk;
%     end
    
    % Crop vRecTemp
    %vRecTemp = vRecTemp(newStart:newEnd, :);% C By Hari
    
    % Final deltaFoF and spks
    tempDeltaFoF = deltaFoF(min(vRecTemp(:,12)):max(vRecTemp(:,12)),:);
    tempSpks = spks(min(vRecTemp(:,12)):max(vRecTemp(:,12)),:);
    fDeltaFoF = [fDeltaFoF; tempDeltaFoF];
    fSpks = [fSpks; tempSpks];
    
    % Add new first and last Ca frame of the trial to lastFrame
    lastFrame(t+1,2) = min(vRecTemp(:,12));
    lastFrame(t+1,3) = max(vRecTemp(:,12));
    
    % Final speed
    tempSpeed = [];
    %%added by Hari from 289-292
    vRecTemp12=vRecTemp(:,12);
    vRecTemp12(isnan(vRecTemp(:,12)))=[];
    vRecTemp2=vRecTemp(:,2);
    vRecTemp2=vRecTemp2(1:length(vRecTemp12));

   % tempSpeed = accumarray(vRecTemp(:,12), vRecTemp(:,2), [],
   % @nanmean);%%original
     tempSpeed = accumarray(vRecTemp12, vRecTemp2, [], @nanmean);
    tempSpeed = tempSpeed(min(vRecTemp(:,12)):end);
    fSpeed = [fSpeed; tempSpeed];
    
    % Final whisk 1
    %%C by Hari
%     tempWhisk1Orig = [];
%     tempWhisk1Orig = accumarray(vRecTemp(:,12), vRecTemp(:,6), [], @nanmean);
%     tempWhisk1Orig = tempWhisk1Orig(min(vRecTemp(:,12)):end);
%     fWhisk1Orig = [fWhisk1Orig; tempWhisk1Orig];
    
    % Final whiskClean 1 and 2
    %c by Hari
%     tempWhisk1 = [];
%     tempWhisk1 = accumarray(vRecTemp(:,12), vRecTemp(:,8), [], @nanmean);
%     tempWhisk1 = tempWhisk1(min(vRecTemp(:,12)):end);
%     fWhisk1 = [fWhisk1; tempWhisk1];
    %c By Hari
%     tempWhisk2 = [];
%     if exist('whiskTrial2', 'var') == 1
%         tempWhisk2 = accumarray(vRecTemp(:,12), vRecTemp(:,9), [], @nanmean);
%         tempWhisk2 = tempWhisk2(min(vRecTemp(:,12)):end);
%         fWhisk2 = [fWhisk2; tempWhisk2];
%     else
%         fWhisk2 = [];
%     end
        
    % Final stim
    % c By Hari
%     tempStim = [];
%     tempStim = accumarray(vRecTemp(:,12), vRecTemp(:,7), [], @nansum);
%     tempStim = tempStim(min(vRecTemp(:,12)):end);
%     fStim = [fStim; tempStim];
    
    % Final sound
    % c by Hari
%     tempSound = [];
%     tempSound = accumarray(vRecTemp(:,12), vRecTemp(:,10), [], @nansum);
%     tempSound = tempSound(min(vRecTemp(:,12)):end);
%     fSound = [fSound; tempSound];
    
    % Final LED
    % c by Hari
%     tempLED = [];
%     tempLED = accumarray(vRecTemp(:,12), vRecTemp(:,11), [], @nansum);
%     tempLED = tempLED(min(vRecTemp(:,12)):end);
%     fLED = [fLED; tempLED];
    
    % Final LEDCL
    % c by Hari
%     tempLEDCL = [];
%     tempLEDCL = accumarray(vRecTemp(:,12), vRecTemp(:,4), [], @nanmean);
%     tempLEDCL = tempLEDCL(min(vRecTemp(:,12)):end);
%     fLEDCL = [fLEDCL; tempLEDCL];
    
    % For convience, save final data as structure
    fTrial.(fieldname).fDeltaFoF = tempDeltaFoF;
    fTrial.(fieldname).fSpks = tempSpks;
    %fTrial.(fieldname).fWhisk1Orig = tempWhisk1Orig;%c By Hari
    %fTrial.(fieldname).fWhisk1 = tempWhisk1;%c By Hari
    %fTrial.(fieldname).fWhisk2 = tempWhisk2;%c By Hari
    fTrial.(fieldname).fSpeed = tempSpeed;
    %fTrial.(fieldname).fStim = tempStim;%c By Hari
    %fTrial.(fieldname).fSound = tempSound;%c By Hari
    %fTrial.(fieldname).fLED = tempLED;%c By Hari
    %fTrial.(fieldname).fLEDCL = tempLEDCL;%c By Hari
end

% --- Save vars ---
%save([currExpVars.folder '\' currExpVars.name], 'caTrialsMask', 'caTrialsMaskLED', 'caTrialsMaskLEDCL', 'lastFrame', 'fDeltaFoF', 'fSpks', 'fWhisk1Orig', 'fWhisk1', 'fWhisk2', 'fSpeed', 'fStim', 'fSound', 'fLED', 'fLEDCL', 'fTrial', '-append') %, '-v7.3')%original
save([currExpVars.folder '\' currExpVars.name], 'caTrialsMask', 'caTrialsMaskLED', 'lastFrame', 'fDeltaFoF', 'fSpks', 'fWhisk1Orig', 'fWhisk1', 'fWhisk2', 'fSpeed', 'fStim', 'fSound', 'fLED', 'fLEDCL', 'fTrial', '-append') %, '-v7.3')%modefied by Hari
% --- Produce raster plot ---
h = figure('units','normalized','outerposition',[0 0 1 1]);

subplot(10,1,[1 6])
climits = [0 10]; % set colorbar limits to be appropriate
iscell = iscell; ce = find(iscell(:,1) == 1); fDeltaFoF1 = fDeltaFoF(:,ce);
imagesc(fDeltaFoF1',climits)
c = colorbar;
pos = get(c,'Position');
c.Position = [pos(1)+.04, pos(2), 0.01, pos(4)]; % move colorbar so graphs are aligned
c.Label.String = 'dF Z-score)';
colormap bone
title('\fontsize{11}dF')
ylabel('Cell #')
xticks(3600);
xticklabels('2 min')
%comment by Hari 384-393
% subplot(10,1,7)
% plot(fWhisk1)
% hold on
% plot(fWhisk2)
% xlim([1 size(fWhisk1,1)]);
% ylim([0 max(fWhisk1)+0.3]);
% title('\fontsize{11}Whisking')
% box off
% xticks(3600);
% xticklabels('\fontsize{10}2 min')

subplot(10,1,8)
plot(fSpeed, 'Color', [.85 .33 .1], 'LineWidth', 1)
xlim([1 size(fSpeed,1)]);
ylim([0 max(fSpeed)+0.3]);
title('\fontsize{11}Locomotion')
ylabel('\fontsize{10}Speed')
box off
xticks(3600);
xticklabels('\fontsize{10}2 min')
%c by Hari 405-425
% subplot(10,1,9)
% plot(fStim/max(fStim), 'Color', [1 .07 .65], 'LineWidth', 1)
% hold on
% plot(fSound/max(fSound), 'Color', [0.3 .075 .93], 'LineWidth', 1)
% xlim([1 size(fStim,1)]);
% ylim([0 1+0.3]);
% title('\fontsize{11}Whisker stim (pink) OR Stimulator Sound (purple)')
% box off
% xticks(3600);
% xticklabels('\fontsize{10}2 min')
% 
% subplot(10,1,10)
% plot(fLED/max(fLED), 'Color', [1 .07 .65], 'LineWidth', 1)
% hold on
% plot(fLEDCL, 'Color', [0.3 .075 .93], 'LineWidth', 1)
% xlim([1 size(fLED,1)]);
% ylim([0 1+0.3]);
% title('\fontsize{11}LED open-Loop (pink) OR LED closed-loop (purple)')
% box off
% xticks(3600);
% xticklabels('\fontsize{10}2 min')

[~, figName, ~] = fileparts([currExpVars.folder '\' currExpVars.name]);
saveas(h, [currExpVars.folder '\' figName '_Raster'])
end