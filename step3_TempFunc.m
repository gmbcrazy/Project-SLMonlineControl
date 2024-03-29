function [fSpeed,timeStampCa_Plane]=PV_SpeedExtract(confSet)


Temp=dir([confSet.save_path0 'TSeries*001*'])
for i=1:length(Temp)
    if a(i).isdir==1
       FilePath=[Temp(i).folder '\' Temp(i).name '\'];
       TrialName=Temp(i).name;
       break
    end
end
info = xml2struct([FilePath  TrialName '.xml']);

numPlanes=length(confSet.ETL);

            timeStampCaTemp = [];
            if numPlanes == 1
                testError = size(info.PVScan.Sequence,2);
                if testError == 1
                    for iFrame = 1:numel(info.PVScan.Sequence.Frame) % Number of frames
                        timeStampCaTemp(iFrame,1) = str2num(info.PVScan.Sequence.Frame{1, iFrame}.Attributes.relativeTime);
                        timeStampCaTemp(iFrame,2) = str2num(info.PVScan.Sequence.Frame{1, iFrame}.Attributes.absoluteTime); % in sec
                        % lastFrame(1,:) = length(timeStampCa);
                    end
                end
                if testError >1
                    for iFrame = 1:numel(info.PVScan.Sequence{1,1}.Frame) % Number of frames
                        timeStampCaTemp(iFrame,1) = str2num(info.PVScan.Sequence{1,1}.Frame{1, iFrame}.Attributes.relativeTime);
                        timeStampCaTemp(iFrame,2) = str2num(info.PVScan.Sequence{1,1}.Frame{1, iFrame}.Attributes.absoluteTime); % in sec
                        % lastFrame(1,:) = length(timeStampCa);
                    end
                end
            else
                for iFrame = 1:numel(info.PVScan.Sequence) % Number of cycles, frame "selectedPlane" within each cycle
                    for selectedPlane=1:numPlanes;
                    timeStampCaTemp(iFrame,1,selectedPlane) = str2num(info.PVScan.Sequence{1, iFrame}.Frame{1, selectedPlane}.Attributes.relativeTime);
                    timeStampCaTemp(iFrame,2,selectedPlane) = str2num(info.PVScan.Sequence{1, iFrame}.Frame{1, selectedPlane}.Attributes.absoluteTime); % in sec
                       if size(info.PVScan.Sequence{1, iFrame}.Frame{1, selectedPlane}.PVStateShard.PVStateValue,2)>1
                          timeStampCaTemp(iFrame,3,selectedPlane) = str2num(info.PVScan.Sequence{1, iFrame}.Frame{1, selectedPlane}.PVStateShard.PVStateValue{1, 1}.Attributes.value);
                       else
                          timeStampCaTemp(iFrame,3,selectedPlane) = str2num(info.PVScan.Sequence{1, iFrame}.Frame{1, selectedPlane}.PVStateShard.PVStateValue.Attributes.value); % I found this issue in the metadata in that the fields are not consistent (in some experiments) -> why?
                       end
                    end
                end
            end

 lastFrame=numel(info.PVScan.Sequence);

 timeStampCa_Plane=squeeze(timeStampCaTemp(:,1,:));


        % --- Read voltage recordings metadata ---
 vFile = dir([FilePath  TrialName '*VoltageRecording*.csv']);
 vRecRaw = csvread([vFile.folder '\' vFile.name],1,0);
 vOutFile = dir([FilePath  TrialName '*VoltageOutput*.xml']);
 vOutTemp = xml2struct([vOutFile.folder '\' vOutFile.name]);
 vOut=vOutTemp;

 vRecXML = dir([FilePath  TrialName '*VoltageRecording*.xml']);
 vRecInfo = xml2struct([vRecXML.folder '\' vRecXML.name]);


 clear whiskTrial

for iPlane=1:numPlanes 
% ---
    vRec=vRecRaw;

    vRecTemp = vRec;
    
    % Delete column 2
    % if deleteColumn2 == 1
    %     vRecTemp = vRecTemp(:,[1 3 4]);
    % end
    
    % Add more working columns
    if size(vRec,2) == 3 % Time, Speed, Cam
        vRec(:,4:12) = NaN; vRec(:,[7 10 11]) = 0;
    end
    if size(vRecTemp,2) == 4 % Time, Speed, Cam, LED
        vRecTemp(:,5:12) = NaN; vRecTemp(:,[7 10 11]) = 0;
    end
    % if exist('BehavStruc')
    %    BehavNum=size(BehavStruc(t).BehData,2);
    %    vRecTemp(:,13:12+BehavNum)=NaN;
    % end
    
    % % Clear column 4 if LEDCL was not on
    % %(sometimes I forgot to disconnect the feedback cable)
    % if caTrialsMaskLEDCL(t) == 0
    %     vRecTemp(:,4) = 0;
    % end
    
    
    sampFreq = 1/diff(vRecTemp(1:2,1))*1000;   %%voltage recording sampling time is micro-seconds
    if sampFreq==str2num(vRecInfo.VRecSessionEntry.Experiment.Rate.Text)

    else
       disp('sampleRate of Voltage Recording does not match;');
    end


    % Cam rising edges matched to nFramesWhisk per trial
    if sum(vRecTemp(:,3)<1)>sum(vRecTemp(:,3)>2)
       IsPositiveTTL=1;
       Temp=-vRecTemp(:,3)+5;
       fallingEdges = LocalMinima(diff(Temp),1,-1)+1; %%TTL pulse, looks like it is negative TTL pulse.
    else
       IsPositiveTTL=0;
       fallingEdges = LocalMinima(diff(vRecTemp(:,3)),1,-1)+1; %%TTL pulse, looks like it is negative TTL pulse.
    end

    % plot(vRecTemp(1:10000,3));hold on;
    % plot(fallingEdges(1:50),vRecTemp(fallingEdges(1:50),3),'r.')
    
    % Indt=find(WhiskerFileID==vRecFileID(t));
    % if ~isempty(Indt)   %% There is possible that some recording are without whisker files recorded.
    if exist('whiskTrial1')
       fW = numel(find(whiskTrial1>0))+1; % Total number of whisk frames actually recorded per trial; +1 because frame 1 is NaN
       fW=min([fW,numel(fallingEdges)]);

    else
       fW=numel(fallingEdges);
    end
    vRecTemp(fallingEdges(1:fW),5) = 1:fW;


    if exist('BehavStruc')
    vRecTemp(fallingEdges(1:fW),13:12+BehavNum)=BehavStruc(Indt).BehData(1:fW,:);
    end
    

    % Add calcium time stamps and frame number
    timeStampCaTrial =  timeStampCa_Plane(:,iPlane)*1000;    %%voltage recording sampling time is micro-seconds
    CaSampleTime=median(diff(timeStampCaTrial));   %%micro-seconds
    % timeStampCaTrial(length(timeStampCaTrial)+1) = timeStampCa(lastFrame(t+1))*1000+mean(diff(timeStampCaTrial(:,1)));
    
    addTime = timeStampCaTrial(1:end);
    addTime(:,2:size(vRecTemp,2)) = NaN;
    vRecTemp = sortrows([addTime; vRecTemp]);  %% Merge Ca signal time and Voltage time togoether
    clear addTime
    
    % fCa = lastFrame(t+1)-(lastFrame(t)+1)+1;  %% # of frames of the trial
    frame = 1;                   %% starting frame Index of the trial
    [~,idx] = intersect(vRecTemp(:,1),timeStampCaTrial(:),'stable');

    if length(idx)<lastFrame+1
       idx(end+1)=idx(end)+length(find(vRecTemp(idx(end)+1:end,1)>=timeStampCaTrial(end)&vRecTemp(idx(end)+1:end,1)<=timeStampCaTrial(end)+CaSampleTime));
    end
    % if fR > 25 NEED TO WORK ON THE MULTIPLANE PART
    for n = 1:lastFrame
            vRecTemp(idx(n):idx(n+1)-1,12) = frame;
            frame = frame + 1;
    end

    % vRecTemp(vRecTemp(:,12)==0,:)=nan;


    % else % NEED TO WORK ON THE MULTIPLANE PART
    %    for n = 1:length(timeStampCaTrial)
    %        vRecNew(timeStampCaTrial(n,1)+1:(timeStampCaTrial(n,1)+1+timeStampCaTrial(n,3)),6) = frame;
    %        frame = frame + 1;
    %    end
    % end
    
    % Establish usable trial time
    startCa = min(find(vRecTemp(:,12) == 1));
    endCa = max(find(vRecTemp(:,12) == lastFrame));
    % startWhisk = find(vRecTemp(:,5) == 1);
    % endWhisk = find(vRecTemp(:,5) == max(vRecTemp(:,5)));
    % 
    % if startWhisk > startCa % means that Whisk aquisition starts after Ca aquisition -> thus is limiting factor; typical stituation except in volume imaging
    %     newStart = startWhisk;
    % else
    %     newStart = startCa;
    % end
    % 
    % if endWhisk > endCa % means that Ca aquisition ended before Whisk aquisition -> thus is limiting factor; typical stituation except in some first experiments in which I recorded only 37500 Whisk frames
    %     % if fR > 25 NEED TO WORK ON THE MULTIPLANE PART
    %         newEnd = endCa;
    %     % else
    %     %    newEnd = endCa + timeStampCaTrial(end,3); NEED TO WORK ON THE MULTIPLANE PART
    %     % end
    % else
    %     newEnd = endWhisk;
    % end
    % if isempty(Indt)  %%in case of no whisker file recorded
    %    newStart = startCa;
    %    newEnd = endCa;
    % end
     newStart(iPlane) = startCa;
     newEnd(iPlane) = endCa;



    % Crop vRecTemp
    vRecTemp = vRecTemp(newStart:newEnd, :);
    
    % Final deltaFoF and spks
    % tempDeltaFoF = deltaFoF(min(vRecTemp(:,12)):max(vRecTemp(:,12)),:);
    % tempSpks = CaData.spks(min(vRecTemp(:,12)):max(vRecTemp(:,12)),:);
    % fDeltaFoF = [fDeltaFoF; tempDeltaFoF];
    % fSpks = CaData.spiks;
    
    % Add new first and last Ca frame of the trial to lastFrame
    % lastFrame(t+1,2) = min(vRecTemp(:,12));
    % lastFrame(t+1,3) = max(vRecTemp(:,12));
    % fSpeed=[];
    % Final speed
    vRecTemp(vRecTemp(:,12)==0,:)=[];
    tempSpeed = [];
    tempSpeed = accumarray(vRecTemp(:,12), vRecTemp(:,2), [], @nanmean);
    tempSpeed = tempSpeed(min(vRecTemp(:,12)):end);
    % fSpeed = [fSpeed; tempSpeed];
    fSpeed(:,iPlane)=tempSpeed;
end
timeStampCa_Plane;
