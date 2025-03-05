PVName='TSeries-03142024-0925-003'       
seqfile='E:\LuRecording\Camera 1\03142024\whisk_03-14-24_11-27-41_01.seq';
nROW=1;  %number of Region of Whisker 
rawFileFolder='E:\LuRecording\03142024_MouseMei03_VGlutAi203SLM\';

rawDataFolder=[rawFileFolder PVName '\'];
xmlFile=[rawDataFolder PVName '.xml']
vFile = dir([rawDataFolder '*VoltageRecording*.csv']);
vFile = [vFile.folder '\' vFile.name];
vOutFile = dir([rawDataFolder '*VoltageOutput*.xml']);
vOutFile=[vOutFile.folder '\' vOutFile.name];

vRec.(fieldname) = csvread([pathGet '\' caTrialCurr.name '\' vFile.name],1,0);
        
        % --- Read voltage output metadata ---
        % --- Use this data to confirm master excel table ---
        vOutFile = dir([pathGet '\' caTrialCurr.name '\*VoltageOutput*.xml']);



[infoTemp, framesWhiskTemp] = Norpix2MATLABAna(seqfile, 100, []);
        % [infoTemp, framesWhiskTemp] = Norpix2MATLABAna([whiskTrialCurr.folder '\' whiskTrialCurr.name], 500, []);
        % implay(framesWhiskTemp,10);
if isfield(infoTemp,'CamFrameRate')
   CamFrameRate=infoTemp.CamFrameRate;
elseif isfield(infoTemp,'FrameRate')
    CamFrameRate=infoTemp.FrameRate;
else
   disp('Check the infoTemp variable, Frame Rate Infomartion missing')
end
framesWhisk = framesWhiskTemp;
totalFrames = infoTemp.AllocatedFrames;
totalFiles=1;


    clear framesWhiskTemp
    framesWhisk = double(framesWhisk);
    posi = [];
   
    if ~isempty(framesWhisk)

    figure; imshow(std(framesWhisk,1,3), [0 5])
    % figure; imagesc(mean(framesWhisk,3))

    % figure; imshow(mean(framesWhisk,3), [0 20])

    for r = 1:nROW
        h = drawrectangle; % should probably use polygon
        pause
        disp 'Press enter when ready'
        posi(r, :) = h.Position;
    end
    close all
    end



[~, framesWhisk] = Norpix2MATLABAna(seqfile, [], []);
framesWhisk=double(framesWhisk);   %%Lu added, change unit8 to double        

tic
for r = 1:size(posi,1)
            % Crop images using the manually set ROI
    r1=max(round(posi(r,2)),1);
    r2=min(round(posi(r,2))+round(posi(r,4)),size(framesWhisk,1));

    r3=max(round(posi(r,1)),1);
    r4=min(round(posi(r,1))+round(posi(r,3)),size(framesWhisk,2));
           
    framesWhiskTemp = framesWhisk(r1:r2, r3:r4, :);


    ft=diff(framesWhiskTemp,1,3);     %%%This way is faster
    ft2 = sqrt(squeeze(sum(sum(ft.^2))));
    ft2=[nan;ft2(:)];
    whiskTrial = ft2;
    % whiskStruc.(fieldname)(:,:,r) = ft2;

end
toc

caTrials=XMLread_SLM(rawFileFolder,3);


infoTrial.(fieldname) = xml2struct([pathGet '\' caTrialCurr.name '\' caTrialCurr.name '.xml']);
        % --- Read voltage recordings metadata ---
vFile = dir([pathGet '\' caTrialCurr.name '\*VoltageRecording*.csv']);
vRec.(fieldname) = csvread([pathGet '\' caTrialCurr.name '\' vFile.name],1,0);
        
        % --- Read voltage output metadata ---
        % --- Use this data to confirm master excel table ---
        infoTrial = xml2struct(xmlFile);
    % We assume that numPlanes, num_x and num_y are constant across the
    % experiment (these parameters MUST be constant across trials)
    numPlanes = 1;
    if size(infoTrial.PVScan.Sequence,2) > 1
        numPlanes = size(infoTrial.PVScan.Sequence{1, 1}.Frame,2);
    end

    timeStampCaTemp(:,1)=caTrials.FrameTS.relativeTime
    timeStampCaTemp(:,2)=caTrials.FrameTS.absoluteTime
    PlaneICaTemp(:,1)=caTrials.FrameTS.relativeTime


    for i=1:length(infoTrial.PVScan.PVStateShard.PVStateValue)

        if isfield(infoTrial.PVScan.PVStateShard.PVStateValue{i}.Attributes,'key')
          i
            infoTrial.PVScan.PVStateShard.PVStateValue{i}.Attributes.key
           if strmatch(infoTrial.PVScan.PVStateShard.PVStateValue{i}.Attributes.key,'linesPerFrame')
              num_x=str2num(infoTrial.PVScan.PVStateShard.PVStateValue{i}.Attributes.value)
           elseif  strmatch(infoTrial.PVScan.PVStateShard.PVStateValue{i}.Attributes.key,'pixelsPerLine')
              num_y=str2num(infoTrial.PVScan.PVStateShard.PVStateValue{i}.Attributes.value)
           else

           end
        end
    end

    vOutTemp = xml2struct(vOutFile);

        % caTrials.LED = vOutTemp.Experiment.Waveform{1,2}.Enabled.Text;
        % caTrials.Evoked = vOutTemp.Experiment.Waveform{1,3}.Enabled.Text;
        % 
        % caTrials.LED_Delay = vOutTemp.Experiment.Waveform{1,2}.WaveformComponent_PulseTrain.FirstPulseDelay.Text;
        % caTrials.LED_Num = vOutTemp.Experiment.Waveform{1,2}.WaveformComponent_PulseTrain.PulseCount.Text;
        % caTrials.LED_Width = vOutTemp.Experiment.Waveform{1,2}.WaveformComponent_PulseTrain.PulseWidth.Text;
        % caTrials.LED_ISI = vOutTemp.Experiment.Waveform{1,2}.WaveformComponent_PulseTrain.PulseSpacing.Text;
        % 
        % caTrials.Evoked_Delay = vOutTemp.Experiment.Waveform{1,3}.WaveformComponent_PulseTrain.FirstPulseDelay.Text;
        % caTrials.Evoked_Num = vOutTemp.Experiment.Waveform{1,3}.WaveformComponent_PulseTrain.PulseCount.Text;
        % caTrials.Evoked_Width = vOutTemp.Experiment.Waveform{1,3}.WaveformComponent_PulseTrain.PulseWidth.Text;
        % caTrials.Evoked_ISI = vOutTemp.Experiment.Waveform{1,3}.WaveformComponent_PulseTrain.PulseSpacing.Text;


PVstateShard=infoTrial.PVScan.PVStateShard.PVStateValue
for i=1:length(PVstateShard)
    if iscell(PVstateShard)
    Temp1=PVstateShard{i};
    elseif isstruct(PVstateShard)
    Temp1=PVstateShard(i);
    else
    end
    
    TempTable1=Temp1.Attributes;

    if isfield(Temp1,'VoltageRecording')
       % struct2table(Temp1.VoltageRecording.Attributes);
       VolFileInfo(i)=Temp1.VoltageRecording.Attributes;
    end
    if isfield(Temp1,'MarkPoints')
       % struct2table(Temp1.VoltageRecording.Attributes);
       % TempTable1=struct2table(Temp1.MarkPoints.Attributes);
       TempTable1=Temp1.MarkPoints.Attributes;
       if ~isempty(TempTable1.filename)
         
       end

       MPrecord(i)=1;
       if  exist('MPFileInfo')
           % MPFileInfo= concatenateStructs(MPFileInfo,TempTable1);
           % MPFileInfo= [MPFileInfo,TempTable1];
            % MPFileInfo = catstruct(MPFileInfo,TempTable1)
            MPFileInfo = concatenateStructs(MPFileInfo,TempTable1);
       else
           MPFileInfo=TempTable1;
       end
    else
        MPrecord(i)=0;
        % MPFileInfo=[];
    end
    if isfield(Temp1,'Frame')
       clear TempframeTS
       numFrames(i)=length(Temp1.Frame);
       for iS=1:length(Temp1.Frame)
           TempframeTS(iS)=Temp1.Frame{iS}.Attributes;           
       end
       if exist('FrameTS')
          FrameTS=[FrameTS;struct2table(TempframeTS)];
       else
          FrameTS=struct2table(TempframeTS);
       end
    else
       numFrames(i)=0;
    end

end

