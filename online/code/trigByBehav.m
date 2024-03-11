%%
% --- 09/23/2022 --- Just reprocessed a few chronic experiments to save a
% much larger triggered block for the manuscript figures, and also save
% trigggered raw/Z2 traces
% --- 08/19/2022 --- Minor spelling
% --- Ana revision July --- In progress; just added a new measurement of
% activity for spontAll and spontLedCl (commented)
% --- Ana revision 06/25/2022 --- around line 509 added if ~isempty
% --- Ana revision 06/16/2022 --- removed the open data on each cycle
% --- Ana revision 05/30/2022 --- no change
% --- Ana revision 05/23/2022 --- major
% --- Ana revision 05/12/2022 --- xlsread replacement
% --- Ana revision 03/21/2022 --- major
% --- Ana revision 02/07/2022 --- Corrected minor issue lines 484-485
% (added -1)
% --- Ana 12302021 --- Initial

function [statsBehSum, statsCaSum, statsCaStarter, rLR] = trigByBehav(currExpVars, baseSave, activity, fR, goodCells, fTrial, fTrialBehav1, fTrialBehav2, ~, groupName, group, paramTrig, starterID)

nCells = numel(goodCells);

%% --- Set parameters for triggering ---
if isempty(paramTrig.win1)
    paramTrig.win1 = round(fR*2);
end
paramTrig.win1New = round(fR/2);
if isempty(paramTrig.win2)
    paramTrig.win2 = round(fR*10); % round(fR);
end
if isempty(paramTrig.win3)
    paramTrig.win3 = round(fR*20); % round(fR*6);
end

%% Extract firing rates across behavioral states
% Should have binarized behavioral state before, during extraction of
% whisking and locomotion onsets-offset detection
movBinary    = [];
groupSub     = [];
groupSubZ    = [];
groupSpks    = [];
groupSpksZ   = [];
groupSpksZ2  = [];
groupSpksMax = [];
for t = 1:size(group,1)
    
    field = group{t};
    
    temp = zeros(size(fTrial.(field).fDeltaFoF,1),1);
    
    for event = 1:numel(fTrialBehav1.(field).onsets1)
        onOnly = fTrialBehav1.(field).onsets1(event)-1;
        if onOnly == 0
            onOnly = 1;
        end
        offOnly = fTrialBehav1.(field).offsets1(event)+1;
        if offOnly >= size(fTrial.(field).fDeltaFoF,1)
            offOnly = size(fTrial.(field).fDeltaFoF,1);
        end
        temp(onOnly:offOnly,1) = 1;
    end
    
    for event = 1:numel(fTrialBehav1.(field).onsets2)
        onOnly = fTrialBehav1.(field).onsets2(event)-1;
        if onOnly == 0
            onOnly = 1;
        end
        offOnly = fTrialBehav1.(field).offsets2(event)+1;
        if offOnly>= size(fTrial.(field).fDeltaFoF,1)
            offOnly = size(fTrial.(field).fDeltaFoF,1);
        end
        if fTrialBehav1.(field).runInd(event) == 1
            temp(onOnly:offOnly,1) = 2;
        end
    end
    movBinary    = [movBinary; temp];
    groupSub     = [groupSub; fTrial.(field).fSub];
    groupSubZ    = [groupSubZ; fTrial.(field).fSubZ];
    groupSpks    = [groupSpks; fTrial.(field).fSpks];
    groupSpksZ   = [groupSpksZ; fTrial.(field).fSpksZ];
    groupSpksZ2  = [groupSpksZ2; fTrial.(field).fSpksZ2];
    groupSpksMax = [groupSpksMax; fTrial.(field).fSpksMax];
end
clear onOnly offOnly

% Establish activity choice (this is all still at prototype level)
rates = [];
if activity == 0
    rates(:,1) = 1:nCells;
    % Firing rate Whole
    rates(:,2)  = mean(groupSpks);
    % Firing rate Q
    rates(:,3)  = mean(groupSpks(movBinary==0,:));
    % Firing rate W
    rates(:,4)  = mean(groupSpks(movBinary==1,:));
    % Firing rate WL
    rates(:,5)  = mean(groupSpks(movBinary==2,:));
    
    % Firing rate Whole
    rates(:,6)  = mean(groupSpksZ);
    % Firing rate Q
    rates(:,7)  = mean(groupSpksZ(movBinary==0,:));
    % Firing rate W
    rates(:,8)  = mean(groupSpksZ(movBinary==1,:));
    % Firing rate WL
    rates(:,9)  = mean(groupSpksZ(movBinary==2,:));
    
    % Firing rate Whole
    rates(:,10) = mean(groupSpksZ2);
    % Firing rate Q
    rates(:,11) = mean(groupSpksZ2(movBinary==0,:));
    % Firing rate W
    rates(:,12) = mean(groupSpksZ2(movBinary==1,:));
    % Firing rate WL
    rates(:,13) = mean(groupSpksZ2(movBinary==2,:));
    
    % Firing rate Whole
    rates(:,14) = mean(groupSpksMax);
    % Firing rate Q
    rates(:,15) = mean(groupSpksMax(movBinary==0,:));
    % Firing rate W
    rates(:,16) = mean(groupSpksMax(movBinary==1,:));
    % Firing rate WL
    rates(:,17) = mean(groupSpksMax(movBinary==2,:));
    clear groupSpks groupSpksZ groupSpksZ2 groupSpksMax t field temp event on off
elseif activity == 1
    rates(:,1) = 1:nCells;
    % Firing rate Whole
    rates(:,2)  = mean(groupSub); % Should probably use integral for DeltaFoF and DeltaFoFZ
    % Firing rate Q
    rates(:,3)  = mean(groupSub(movBinary==0,:));
    % Firing rate W
    rates(:,4)  = mean(groupSub(movBinary==1,:));
    % Firing rate WL
    rates(:,5)  = mean(groupSub(movBinary==2,:));
    
    % Firing rate Whole
    rates(:,6)  = mean(groupSubZ);
    % Firing rate Q
    rates(:,7)  = mean(groupSubZ(movBinary==0,:));
    % Firing rate W
    rates(:,8)  = mean(groupSubZ(movBinary==1,:));
    % Firing rate WL
    rates(:,9)  = mean(groupSubZ(movBinary==2,:));
    
    clear groupDeltaFoF groupdeltaFoFZ t field temp event on off
end

%% Calculate R values
% Concatenate trials
corrCa = [];
corrWhisk1 = [];
corrWhisk2 = [];
corrSpeed = [];
for t = 1:size(group,1)
    field = group{t};
    corrCa = [corrCa; fTrial.(field).choice1];
    corrWhisk1 = [corrWhisk1; fTrialBehav1.(field).fWhisk1Norm];
    if ~isempty(fTrialBehav2)
        corrWhisk2 = [corrWhisk2; fTrialBehav2.(field).fWhisk2Norm];
    end
    corrSpeed = [corrSpeed; fTrial.(field).fSpeed];
end
clear t field

% Remove NaN elements (could do it directly using corr)
corrCa(isnan(corrWhisk1),:) = [];
corrSpeed(isnan(corrWhisk1)) = [];
if ~isempty(fTrialBehav2)
    corrWhisk2(isnan(corrWhisk1)) = [];
end
corrWhisk1(isnan(corrWhisk1)) = [];

% Remove high frequency oscillations from whisk1 and 2
corrWhisk1Upper = smoothdata(corrWhisk1, 'movmean', 10);
% [corrWhisk1Upper,~] = envelope(corrWhisk1);
if ~isempty(fTrialBehav2)
    [corrWhisk2Upper,~] = smoothdata(corrWhisk2, 'movmean', 10);
end

% If avaliable, correlated left and right whisker motion
rLR = [];
if ~isempty(fTrialBehav2)
    [rLR(:,1), ~] = corr(corrWhisk1, corrWhisk2);
    [rLR(:,2), ~] = corr(corrWhisk1Upper, corrWhisk2Upper);
else
    rLR(1,1:2) = NaN;
end

% No binning (for comparison)
RWhiskSpeed = [];
RWhiskSpeed(:,1) = 1:nCells;

% [RWhiskSpeed(:,2),~] = corr(corrCa, corrWhisk1); % no need to use Z-scoring for corr
% [RWhiskSpeed(:,3),~] = corr(corrCa, corrWhisk1, 'Type', 'Spearman');
[RWhiskSpeed(:,2),~] = corr(corrCa, corrWhisk1Upper);
% [RWhiskSpeed(:,3),~] = corr(corrCa, corrWhisk1Upper, 'Type', 'Spearman');
[RWhiskSpeed(:,3),~] = corr(corrCa, corrSpeed);
% [RWhiskSpeed(:,5),~] = corr(corrCa, corrSpeed, 'Type', 'Spearman');

% Bin fSpks and behavior (does not make a lot of sense to bin DeltaFoF, but
% just to follow the same routine and make everything consistent)
counter = 3;
b = [6 2 1];
for bb = 1:size(b,2)
    
    x_points = round(fR/b(bb));
    edges = 0:x_points:size(corrCa,1);
    
    newCorrCa = [];
    newCorrWhisk1 = [];
    newCorrSpeed = [];
    for e = 1:numel(edges)-1
        newCorrCa = [newCorrCa; sum(corrCa(edges(e)+1:edges(e+1),:))];
        newCorrWhisk1 = [newCorrWhisk1; sum(corrWhisk1(edges(e)+1:edges(e+1),1))];
        newCorrSpeed = [newCorrSpeed; sum(corrSpeed(edges(e)+1:edges(e+1),1))];
    end
    
    [RWhiskSpeed(:,counter+1),~] = corr(newCorrCa, newCorrWhisk1);
    % [RWhiskSpeed(:,counter+2),~] = corr(newCorrCa, newCorrWhisk1, 'Type', 'Spearman');
    [RWhiskSpeed(:,counter+2),~] = corr(newCorrCa, newCorrSpeed);
    % [RWhiskSpeed(:,counter+4),~] = corr(newCorrCa, newCorrSpeed, 'Type', 'Spearman');
    
    counter = counter+2;
end
clear x_points edges e corrCa corrWhisk corrSpeed newCorrCa newCorrWhisk1 newCorrWhisk1Upper newCorrSpeed b bb counter

%% --- Collect data for plotting, stats and sorting ---
caAvW    = [];
blockWS  = []; blockWE  = [];
caAvWL   = []; caAvDiff = []; caAvSilence = [];
blockWLS = []; blockWLE = [];
statsBeh = [];
eventCounter = 1;
eventCounter0 = 1; eventCounter00 = 1; eventCounter000 = 1;
eventCounter1 = 1; eventCounter11 = 1; eventCounter111 = 1;
for t = 1:size(group,1)
    
    field = group{t};
    
    for event = 1:numel(fTrialBehav1.(field).onsets2)
        
        onOnly = fTrialBehav1.(field).onsets2(event)-1;
        if onOnly == 0
            onOnly = 1;
        end
        offOnly = fTrialBehav1.(field).offsets2(event)+1;
        if offOnly>= size(fTrial.(field).fDeltaFoF,1)
            offOnly = size(fTrial.(field).fDeltaFoF,1);
        end
        
        if fTrialBehav1.(field).runInd(event) == 0
            % durTemp = [0 offOnly-onOnly offOnly-onOnly+1 0 0 0];
            statsBeh(eventCounter, :) = [0 offOnly-onOnly offOnly-onOnly+1 0 0 0];
            
            if onOnly-paramTrig.win1New>0 && offOnly+paramTrig.win1New<=size(fTrial.(field).fDeltaFoF,1)
                % durTemp = [0 offOnly-onOnly offOnly-onOnly+1 0 0 1];
                % statsBeh = [statsBeh; durTemp];
                statsBeh(eventCounter, 6) = 1;
                for cell = 1:nCells
                    % base = mean(fTrial.(field).fDeltaFoF(on-paramTrig.win1:on-1, cell));
                    caAvW(cell,1,eventCounter0)  = mean(fTrial.(field).choice1(onOnly-paramTrig.win1New: onOnly-1, cell));% -base;
                    caAvW(cell,2,eventCounter0)  = mean(fTrial.(field).choice1(onOnly+1: offOnly-1, cell)); % -base;
                    caAvW(cell,3,eventCounter0)  = mean(fTrial.(field).choice1(offOnly+1: offOnly+paramTrig.win1New, cell)); % -base;
                    
                    caAvW(cell,4,eventCounter0)  = mean(fTrial.(field).choice2(onOnly-paramTrig.win1New: onOnly-1, cell));% -base;
                    caAvW(cell,5,eventCounter0)  = mean(fTrial.(field).choice2(onOnly+1: offOnly-1, cell)); % -base;
                    caAvW(cell,6,eventCounter0)  = mean(fTrial.(field).choice2(offOnly+1: offOnly+paramTrig.win1New, cell)); % -base;
                    
                    caAvW(cell,7,eventCounter0)  = caAvW(cell,2,eventCounter0) - caAvW(cell,1,eventCounter0);
                    caAvW(cell,8,eventCounter0)  = caAvW(cell,3,eventCounter0) - caAvW(cell,2,eventCounter0);
                    
                    caAvW(cell,9,eventCounter0)  = caAvW(cell,5,eventCounter0) - caAvW(cell,4,eventCounter0);
                    caAvW(cell,10,eventCounter0) = caAvW(cell,6,eventCounter0) - caAvW(cell,5,eventCounter0);
                    
                    if onOnly-paramTrig.win1>0 && onOnly+paramTrig.win2<=size(fTrial.(field).fDeltaFoF,1)
                        blockWS.Ca(cell,:,eventCounter00) = fTrial.(field).choice2(onOnly-paramTrig.win1: onOnly+paramTrig.win2, cell);
                        blockWS.CaRaw(cell,:,eventCounter00) = fTrial.(field).choice1(onOnly-paramTrig.win1: onOnly+paramTrig.win2, cell);
                        blockWS.CaZ2(cell,:,eventCounter00) = fTrial.(field).choice3(onOnly-paramTrig.win1: onOnly+paramTrig.win2, cell);
                    end
                    if offOnly-paramTrig.win2>0 && offOnly+paramTrig.win1<=size(fTrial.(field).fDeltaFoF,1)
                        blockWE.Ca(cell,:,eventCounter000) = fTrial.(field).choice2(offOnly-paramTrig.win2: offOnly+paramTrig.win1, cell);
                        blockWE.CaRaw(cell,:,eventCounter000) = fTrial.(field).choice1(offOnly-paramTrig.win2: offOnly+paramTrig.win1, cell);
                        blockWE.CaZ2(cell,:,eventCounter000) = fTrial.(field).choice3(offOnly-paramTrig.win2: offOnly+paramTrig.win1, cell);
                    end
                end
                eventCounter0 = eventCounter0+1;
                
                if onOnly-paramTrig.win1>0 && onOnly+paramTrig.win2<=size(fTrial.(field).fDeltaFoF,1)
                    blockWS.Whisk1(:,eventCounter00) = fTrial.(field).fWhisk1(onOnly-paramTrig.win1: onOnly+paramTrig.win2);
                    blockWS.Speed(:,eventCounter00)  = fTrial.(field).fSpeed(onOnly-paramTrig.win1: onOnly+paramTrig.win2);
                    blockWS.LedCl(:,eventCounter00)  = fTrial.(field).fLEDCL(onOnly-paramTrig.win1: onOnly+paramTrig.win2);
                    blockWS.Led(:,eventCounter00)    = fTrial.(field).fLED(onOnly-paramTrig.win1: onOnly+paramTrig.win2);
                    blockWS.Stim(:,eventCounter00)   = fTrial.(field).fStim(onOnly-paramTrig.win1: onOnly+paramTrig.win2);
                    blockWS.Sound(:,eventCounter00)  = fTrial.(field).fSound(onOnly-paramTrig.win1: onOnly+paramTrig.win2);
                    eventCounter00 = eventCounter00+1;
                end
                
                if offOnly-paramTrig.win2>0 && offOnly+paramTrig.win1<=size(fTrial.(field).fDeltaFoF,1)
                    blockWE.Whisk1(:,eventCounter000) = fTrial.(field).fWhisk1(offOnly-paramTrig.win2: offOnly+paramTrig.win1);
                    blockWE.Speed(:,eventCounter000)  = fTrial.(field).fSpeed(offOnly-paramTrig.win2: offOnly+paramTrig.win1);
                    blockWE.LedCl(:,eventCounter000)  = fTrial.(field).fLEDCL(offOnly-paramTrig.win2: offOnly+paramTrig.win1);
                    blockWE.Led(:,eventCounter000)    = fTrial.(field).fLED(offOnly-paramTrig.win2: offOnly+paramTrig.win1);
                    blockWE.Stim(:,eventCounter000)   = fTrial.(field).fStim(offOnly-paramTrig.win2: offOnly+paramTrig.win1);
                    blockWE.Sound(:,eventCounter000)  = fTrial.(field).fSound(offOnly-paramTrig.win2: offOnly+paramTrig.win1);
                    eventCounter000 = eventCounter000+1;
                end
            end
        end
        if fTrialBehav1.(field).runInd(event) == 1
            temp = fTrial.(field).fSpeed(onOnly:offOnly);
            avSpeed = mean(temp(temp>0.025)); % remove points >= 0.025 because speed does not necessarily start with whisk
            maxSpeed = max(temp);
            % durTemp = [1 offOnly-onOnly offOnly-onOnly+1 avSpeed maxSpeed 0];
            statsBeh(eventCounter, :) = [1 offOnly-onOnly offOnly-onOnly+1 avSpeed maxSpeed 0];
            
            if onOnly-paramTrig.win1New>0 && offOnly+paramTrig.win1New<=size(fTrial.(field).fDeltaFoF,1)
                % durTemp = [1 offOnly-onOnly offOnly-onOnly+1 avSpeed maxSpeed 1];
                % statsBeh = [statsBeh; durTemp];
                statsBeh(eventCounter, 6) = 1;
                for cell = 1:nCells
                    % base = mean(fTrial.(field).fDeltaFoF(on-paramTrig.win1:on-1, cell));
                    caAvWL(cell,1,eventCounter1)  = mean(fTrial.(field).choice1(onOnly-paramTrig.win1New: onOnly-1, cell)); % -base;
                    caAvWL(cell,2,eventCounter1)  = mean(fTrial.(field).choice1(onOnly+1: offOnly-1, cell)); % -base;
                    caAvWL(cell,3,eventCounter1)  = mean(fTrial.(field).choice1(offOnly+1: offOnly+paramTrig.win1New, cell)); % -base;
                    
                    caAvWL(cell,4,eventCounter1)  = mean(fTrial.(field).choice2(onOnly-paramTrig.win1New: onOnly-1, cell));% -base;
                    caAvWL(cell,5,eventCounter1)  = mean(fTrial.(field).choice2(onOnly+1: offOnly-1, cell)); % -base;
                    caAvWL(cell,6,eventCounter1)  = mean(fTrial.(field).choice2(offOnly+1: offOnly+paramTrig.win1New, cell)); % -base;
                    
                    caAvWL(cell,7,eventCounter1)  = caAvWL(cell,2,eventCounter1) - caAvWL(cell,1,eventCounter1);
                    caAvWL(cell,8,eventCounter1)  = caAvWL(cell,3,eventCounter1) - caAvWL(cell,2,eventCounter1);
                    
                    caAvWL(cell,9,eventCounter1)  = caAvWL(cell,5,eventCounter1) - caAvWL(cell,4,eventCounter1);
                    caAvWL(cell,10,eventCounter1) = caAvWL(cell,6,eventCounter1) - caAvWL(cell,5,eventCounter1);
                    
%                     temp = fTrial.(field).choice1(onOnly+1: offOnly-1, cell);
%                     ana = find(temp>0);
%                     caAvDiff(eventCounter1,cell) = mean(diff(ana));
%                     ana = find(temp<10);
%                     caAvSilence(eventCounter1,cell) = numel(ana)/numel(temp);

                    if onOnly-paramTrig.win1>0 && onOnly+paramTrig.win3<=size(fTrial.(field).fDeltaFoF,1)
                        blockWLS.Ca(cell,:,eventCounter11) = fTrial.(field).choice2(onOnly-paramTrig.win1: onOnly+paramTrig.win3, cell);
                        blockWLS.CaRaw(cell,:,eventCounter11) = fTrial.(field).choice1(onOnly-paramTrig.win1: onOnly+paramTrig.win3, cell);
                        blockWLS.CaZ2(cell,:,eventCounter11) = fTrial.(field).choice3(onOnly-paramTrig.win1: onOnly+paramTrig.win3, cell);
                    end
                    if offOnly-paramTrig.win3>0 && offOnly+paramTrig.win1<=size(fTrial.(field).fDeltaFoF,1)
                        blockWLE.Ca(cell,:,eventCounter111) = fTrial.(field).choice2(offOnly-paramTrig.win3: offOnly+paramTrig.win1, cell);
                        blockWLE.CaRaw(cell,:,eventCounter111) = fTrial.(field).choice1(offOnly-paramTrig.win3: offOnly+paramTrig.win1, cell);
                        blockWLE.CaZ2(cell,:,eventCounter111) = fTrial.(field).choice3(offOnly-paramTrig.win3: offOnly+paramTrig.win1, cell);
                    end
                end
                eventCounter1 = eventCounter1+1;
                
                if onOnly-paramTrig.win1>0 && onOnly+paramTrig.win3<=size(fTrial.(field).fDeltaFoF,1)
                    blockWLS.Whisk1(:,eventCounter11) = fTrial.(field).fWhisk1(onOnly-paramTrig.win1: onOnly+paramTrig.win3);
                    blockWLS.Speed(:,eventCounter11)  = fTrial.(field).fSpeed(onOnly-paramTrig.win1: onOnly+paramTrig.win3);
                    blockWLS.LedCl(:,eventCounter11)  = fTrial.(field).fLEDCL(onOnly-paramTrig.win1: onOnly+paramTrig.win3);
                    blockWLS.Led(:,eventCounter11)    = fTrial.(field).fLED(onOnly-paramTrig.win1: onOnly+paramTrig.win3);
                    blockWLS.Stim(:,eventCounter11)   = fTrial.(field).fStim(onOnly-paramTrig.win1: onOnly+paramTrig.win3);
                    blockWLS.Sound(:,eventCounter11)  = fTrial.(field).fSound(onOnly-paramTrig.win1: onOnly+paramTrig.win3);
                    eventCounter11 = eventCounter11+1;
                end
                
                if  offOnly-paramTrig.win3>0 && offOnly+paramTrig.win1<=size(fTrial.(field).fDeltaFoF,1)
                    blockWLE.Whisk1(:,eventCounter111) = fTrial.(field).fWhisk1(offOnly-paramTrig.win3: offOnly+paramTrig.win1);
                    blockWLE.Speed(:,eventCounter111)  = fTrial.(field).fSpeed(offOnly-paramTrig.win3: offOnly+paramTrig.win1);
                    blockWLE.LedCl(:,eventCounter111)  = fTrial.(field).fLEDCL(offOnly-paramTrig.win3: offOnly+paramTrig.win1);
                    blockWLE.Led(:,eventCounter111)    = fTrial.(field).fLED(offOnly-paramTrig.win3: offOnly+paramTrig.win1);
                    blockWLE.Stim(:,eventCounter111)   = fTrial.(field).fStim(offOnly-paramTrig.win3: offOnly+paramTrig.win1);
                    blockWLE.Sound(:,eventCounter111)  = fTrial.(field).fSound(offOnly-paramTrig.win3: offOnly+paramTrig.win1);
                    eventCounter111 = eventCounter111+1;
                end
            end
        end
        eventCounter = eventCounter+1;
    end
end
clear eventCounter00 eventCounter000 eventCounter11 eventCounter111 t fieldname event on off cell base temp avSpeed maxSpeed durTemp

%% --- Simple stats ----
statsCa = zeros(nCells,28);
statsCa(:,1) = 1:nCells;
for cell = 1:nCells
    if size(caAvW,3) < 4
        statsCa(cell,2) = NaN;
        statsCa(cell,3) = NaN;
        
        statsCa(cell,6)  = NaN;
        statsCa(cell,7)  = NaN;
        statsCa(cell,8)  = NaN;
        statsCa(cell,9)  = NaN; % Z
        statsCa(cell,10) = NaN; % Z
        statsCa(cell,11) = NaN; % Z
        statsCa(cell,12) = NaN; % ModInd
        statsCa(cell,13) = NaN; % ModInd
        statsCa(cell,14) = NaN; % ModIndZ
        statsCa(cell,15) = NaN; % ModIndZ
    else
        [~, statsCa(cell,2)] = ttest(squeeze(caAvW(cell,2,:)), squeeze(caAvW(cell,1,:)));
        [~, statsCa(cell,3)] = ttest(squeeze(caAvW(cell,3,:)), squeeze(caAvW(cell,2,:)));
        if isnan(statsCa(cell,2)) % I have come across examples in which, with very few behavioral events, activities can be zero across events, generating NaNs; these cells can still be taken as controls
            statsCa (cell,2) = 1;
        end
        if isnan(statsCa(cell,3)) % I have come across examples in which, with very few behavioral events, activities can be zero across events, generating NaNs; these cells can still be taken as controls
            statsCa (cell,3) = 1;
        end
        statsCa(cell,6)  = mean(squeeze(caAvW(cell,1,:)));
        statsCa(cell,7)  = mean(squeeze(caAvW(cell,2,:)));
        statsCa(cell,8)  = mean(squeeze(caAvW(cell,3,:)));
        statsCa(cell,9)  = mean(squeeze(caAvW(cell,4,:))); % Z
        statsCa(cell,10) = mean(squeeze(caAvW(cell,5,:))); % Z
        statsCa(cell,11) = mean(squeeze(caAvW(cell,6,:))); % Z
        statsCa(cell,12) = mean(squeeze(caAvW(cell,7,:))); % ModInd
        statsCa(cell,13) = mean(squeeze(caAvW(cell,8,:))); % ModInd
        statsCa(cell,14) = mean(squeeze(caAvW(cell,9,:))); % ModIndZ
        statsCa(cell,15) = mean(squeeze(caAvW(cell,10,:))); % ModIndZ
    end
    if size(caAvWL,3) < 4
        statsCa(cell,4) = NaN;
        statsCa(cell,5) = NaN;
        
        statsCa(cell,16) = NaN;
        statsCa(cell,17) = NaN;
        statsCa(cell,18) = NaN;
        statsCa(cell,19) = NaN; % Z
        statsCa(cell,20) = NaN; % Z
        statsCa(cell,21) = NaN; % Z
        statsCa(cell,22) = NaN; % ModInd
        statsCa(cell,23) = NaN; % ModInd
        statsCa(cell,24) = NaN; % ModIndZ
        statsCa(cell,25) = NaN; % ModIndZ
    else
        [~, statsCa(cell,4)] = ttest(squeeze(caAvWL(cell,2,:)), squeeze(caAvWL(cell,1,:)));
        [~, statsCa(cell,5)] = ttest(squeeze(caAvWL(cell,3,:)), squeeze(caAvWL(cell,2,:)));
        if isnan(statsCa(cell,4)) % I have come across examples in which, with very few behavioral events, activities can be zero across events, generating NaNs; these cells can still be taken as controls
            statsCa (cell,4) = 1;
        end
        if isnan(statsCa(cell,5)) % I have come across examples in which, with very few behavioral events, activities can be zero across events, generating NaNs; these cells can still be taken as controls
            statsCa (cell,5) = 1;
        end
        statsCa(cell,16)  = mean(squeeze(caAvWL(cell,1,:)));
        statsCa(cell,17)  = mean(squeeze(caAvWL(cell,2,:)));
        statsCa(cell,18)  = mean(squeeze(caAvWL(cell,3,:)));
        statsCa(cell,19)  = mean(squeeze(caAvWL(cell,4,:))); % Z
        statsCa(cell,20)  = mean(squeeze(caAvWL(cell,5,:))); % Z
        statsCa(cell,21)  = mean(squeeze(caAvWL(cell,6,:))); % Z
        statsCa(cell,22)  = mean(squeeze(caAvWL(cell,7,:))); % ModInd
        statsCa(cell,23)  = mean(squeeze(caAvWL(cell,8,:))); % ModInd
        statsCa(cell,24)  = mean(squeeze(caAvWL(cell,9,:))); % ModIndZ
        statsCa(cell,25)  = mean(squeeze(caAvWL(cell,10,:))); % ModIndZ
    end
end
statsCa(:,26) = RWhiskSpeed(:,6);
statsCa(:,27) = RWhiskSpeed(:,7);
% Sort significant cells by Mod ind (just for onsets for the moment)
% Arrangement more convinent for plotting

% unchangW = []; downregW = []; upregW = [];
% unchangWS = statsCa(statsCa(:,2,:)>=0.05,:); unchangWS = sortrows(unchangWS, 12, 'ascend');
% downregWS = statsCa(statsCa(:,2,:)<0.05 & statsCa(:,12,:)<0,:); downregWS = sortrows(downregWS, 12, 'descend');
% upregWS   = statsCa(statsCa(:,2,:)<0.05 & statsCa(:,12,:)>0,:); upregWS = sortrows(upregWS, 12, 'ascend');

% unchangWL = []; downregWL = []; upregWL = [];
% unchangWLS = statsCa(statsCa(:,4,:)>=0.05,:); unchangWLS = sortrows(unchangWLS, 22, 'ascend');
% downregWLS = statsCa(statsCa(:,4,:)<0.05 & statsCa(:,22,:)<0,:); downregWLS = sortrows(downregWLS, 22, 'descend');
% upregWLS   = statsCa(statsCa(:,4,:)<0.05 & statsCa(:,22,:)>0,:); upregWLS = sortrows(upregWLS, 22, 'ascend');

% Arrangement for general stats (rewrite; ugly)
% Wtemp     = [];
% Wtemp     = statsCa(statsCa(:,2,:)<0.05,:);
% Wonly     = Wtemp(Wtemp(:,4,:)>=0.05,:);
% Wboth     = Wtemp(Wtemp(:,4,:)<0.05,:);
%
% Wdown     = Wonly(Wonly(:,12,:)<0,:); % Wdown = sortrows(Wdown, 6, 'descend'); %%%
% Wup       = Wonly(Wonly(:,12,:)>0,:); % Wup = sortrows(Wup, 6, 'ascend'); %%%
%
% Wbothdown = Wboth(Wboth(:,12,:)<0 & Wboth(:,22,:)<0,:); % Wbothdown = sortrows(Wbothdown, 8, 'descend'); %%%
% Wbothup   = Wboth(Wboth(:,12,:)>0 & Wboth(:,22,:)>0,:); % Wbothup = sortrows(Wbothup, 8, 'ascend'); %%%
% Wbothrev1 = Wboth(Wboth(:,12,:)<0 & Wboth(:,22,:)>0,:); %%%
% Wbothrev2 = Wboth(Wboth(:,12,:)>0 & Wboth(:,22,:)<0,:); %%%
%
% Wtemp     = [];
% Wtemp     = statsCa(statsCa(:,2,:)>=0.05 & statsCa(:,4,:)<0.05,:);
%
% WLdown    = Wtemp(Wtemp(:,22,:)<0,:); % WLdown = sortrows(WLdown, 8, 'descend'); %%%
% WLup      = Wtemp(Wtemp(:,22,:)>0,:); % WLup = sortrows(WLup, 8, 'ascend'); %%%

% clear Wtemp Wonly Wboth

% ---
W = []; WCtrl = []; WDown = []; WUp = [];

onAndOff   = statsCa(statsCa(:,2,:)<0.01 & statsCa(:,3,:)<0.01,:);
onOnly     = statsCa(statsCa(:,2,:)<0.01 & statsCa(:,3,:)>=0.01,:);
offOnly    = statsCa(statsCa(:,2,:)>=0.01 & statsCa(:,3,:)<0.01,:);

W.up1      = onAndOff(onAndOff(:,12,:)>0 & onAndOff(:,13,:)<0,:);
W.down1    = onAndOff(onAndOff(:,12,:)<0 & onAndOff(:,13,:)>0,:);

W.strange1 = onAndOff(onAndOff(:,12,:)>0 & onAndOff(:,13,:)>0,:); % It might happen, particularly when using deltaFoF
W.strange2 = onAndOff(onAndOff(:,12,:)<0 & onAndOff(:,13,:)<0,:);

W.up2      = onOnly(onOnly(:,12,:)>0,:);
W.down2    = onOnly(onOnly(:,12,:)<0,:);

W.up3      = offOnly(offOnly(:,13,:)<0,:);
W.down3    = offOnly(offOnly(:,13,:)>0,:);

WCtrl      = statsCa(statsCa(:,2,:)>=0.01 & statsCa(:,3,:)>=0.01,:); WCtrl = sortrows(WCtrl, 14, 'ascend'); %%%
WUp        = [W.up1; W.up2; W.up3; W.strange1]; WUp = sortrows(WUp, 14, 'ascend');
WDown      = [W.down1; W.down2; W.down3; W.strange2]; WDown = sortrows(WDown, 14, 'descend');

% ---
WL = []; WLCtrl = []; WLDown = []; WLUp = [];

onAndOff    = statsCa(statsCa(:,4,:)<0.01 & statsCa(:,5,:)<0.01,:);
onOnly      = statsCa(statsCa(:,4,:)<0.01 & statsCa(:,5,:)>=0.01,:);
offOnly     = statsCa(statsCa(:,4,:)>=0.01 & statsCa(:,5,:)<0.01,:);

WL.up1      = onAndOff(onAndOff(:,22,:)>0 & onAndOff(:,23,:)<0,:);
WL.down1    = onAndOff(onAndOff(:,22,:)<0 & onAndOff(:,23,:)>0,:);

WL.strange1 = onAndOff(onAndOff(:,22,:)>0 & onAndOff(:,23,:)>0,:); % It might happen, particularly when using deltaFoF
WL.strange2 = onAndOff(onAndOff(:,22,:)<0 & onAndOff(:,23,:)<0,:);

WL.up2      = onOnly(onOnly(:,22,:)>0,:);
WL.down2    = onOnly(onOnly(:,22,:)<0,:);

WL.up3      = offOnly(offOnly(:,23,:)<0,:);
WL.down3    = offOnly(offOnly(:,23,:)>0,:);

WLCtrl      = statsCa(statsCa(:,4,:)>=0.01 & statsCa(:,5,:)>=0.01,:); WLCtrl = sortrows(WLCtrl, 24, 'ascend'); %%%
WLUp        = [WL.up1; WL.up2; WL.up3; WL.strange1]; WLUp = sortrows(WLUp, 24, 'ascend');
WLDown      = [WL.down1; WL.down2; WL.down3; WL.strange2]; WLDown = sortrows(WLDown, 24, 'descend');

clear onAndOff onOnly offOnly
% ---
WBothCtrl   = intersect(WCtrl(:,1), WLCtrl(:,1));
WBothUp     = intersect(WUp(:,1), WLUp(:,1));
WBothDown   = intersect(WDown(:,1), WLDown(:,1));
WBothRev1   = intersect(WUp(:,1), WLDown(:,1));
WBothRev2   = intersect(WDown(:,1), WLUp(:,1));

WUpFinal    = intersect(WUp(:,1), WLCtrl(:,1));
WLUpFinal   = intersect(WLUp(:,1), WCtrl(:,1));
WDownFinal  = intersect(WDown(:,1), WLCtrl(:,1));
WLDownFinal = intersect(WLDown(:,1), WCtrl(:,1));

if isnan(statsCa(:,2:3))
    WBothCtrl   = WLCtrl;
    WLUpFinal   = WLUp(:,1);
    WLDownFinal = WLDown(:,1);
end
if isnan(statsCa(:,4:5))
    WBothCtrl   = WCtrl;
    WUpFinal    = WUp(:,1);
    WDownFinal  = WDown(:,1);
end

% statsCa(:,28)              = 0;
statsCa(WDownFinal(:,1),28)  = -1;
statsCa(WUpFinal(:,1),28)    = +1;
statsCa(WLDownFinal(:,1),28) = -2;
statsCa(WLUpFinal(:,1),28)   = +2;
if ~isempty(WBothDown)
    statsCa(WBothDown(:,1),28)   = -3;
end
if ~isempty(WBothUp)
    statsCa(WBothUp(:,1),28)     = +3;
end
if ~isempty(WBothRev1)
    statsCa(WBothRev1(:,1),28) = +4;
end
if ~isempty(WBothRev2)
    statsCa(WBothRev2(:,1),28) = +5;
end
statsCaSum = [];
statsCaSum = [nCells size(WBothCtrl,1) size(WDownFinal,1) size(WUpFinal,1) size(WLDownFinal,1) size(WLUpFinal,1) size(WBothDown,1) size(WBothUp,1) size(WBothRev1,1) size(WBothRev2,1)];
clear cell

%% --- Time normalize ---
blockWLSresized = [];
maxWL = max(statsBeh((statsBeh(:,1)==1),2)); % 4500 (used for reprocessing of some chronic exp Sep 2022 
eventCounter = 1;
for t = 1:size(group,1)
    
    field = group{t};
    
    for event = 1:numel(fTrialBehav1.(field).onsets2)
        
        onOnly = fTrialBehav1.(field).onsets2(event)-1;
        offOnly = fTrialBehav1.(field).offsets2(event)+1;
        
        if onOnly-paramTrig.win2>0 && offOnly+paramTrig.win2<=size(fTrial.(field).fDeltaFoF,1)
            
            if fTrialBehav1.(field).runInd(event) == 1
                x      = 1:1:(offOnly-onOnly+1);
                xi     = linspace(1, (offOnly-onOnly+1), maxWL)';
                
                for cell = 1:nCells
                    middle = interp1(x, fTrial.(field).choice2(onOnly:offOnly, cell), xi);
                    bef    = fTrial.(field).choice2(onOnly-paramTrig.win2:onOnly-1, cell);
                    aft    = fTrial.(field).choice2(offOnly+1:offOnly+paramTrig.win2, cell);
                    blockWLSresized.Ca(cell,:,eventCounter) = [bef; middle; aft];
                end
                middle = interp1(x, fTrial.(field).fWhisk1(onOnly:offOnly), xi);
                bef    = fTrial.(field).fWhisk1(onOnly-paramTrig.win2:onOnly-1);
                aft    = fTrial.(field).fWhisk1(offOnly+1:offOnly+paramTrig.win2);
                blockWLSresized.Whisk1(:,eventCounter) = [bef; middle; aft];
                
                middle = interp1(x, fTrial.(field).fSpeed(onOnly:offOnly), xi);
                bef    = fTrial.(field).fSpeed(onOnly-paramTrig.win2:onOnly-1);
                aft    = fTrial.(field).fSpeed(offOnly+1:offOnly+paramTrig.win2);
                blockWLSresized.Speed(:,eventCounter) = [bef; middle; aft];
                
                middle = interp1(x, fTrial.(field).fLED(onOnly:offOnly), xi);
                bef    = fTrial.(field).fLED(onOnly-paramTrig.win2:onOnly-1);
                aft    = fTrial.(field).fLED(offOnly+1:offOnly+paramTrig.win2);
                blockWLSresized.Led(:,eventCounter) = [bef; middle; aft];
                
                middle = interp1(x, fTrial.(field).fLEDCL(onOnly:offOnly), xi);
                bef    = fTrial.(field).fLEDCL(onOnly-paramTrig.win2:onOnly-1);
                aft    = fTrial.(field).fLEDCL(offOnly+1:offOnly+paramTrig.win2);
                blockWLSresized.LedCl(:,eventCounter) = [bef; middle; aft];
                eventCounter = eventCounter+1;
            end
        end
    end
end
clear maxWL eventCounter t field event on off x xi cell middle bef aft

% Find bin with max firing

%% --- Sum up behavior ---
IEI = [];
l = [];
excluded = [];
count = 1;
for t = 1:size(group,1)
    field = group{t};
    l(t) = size(fTrial.(field).fDeltaFoF,1);
    excludedOnset = setdiff(fTrialBehav1.(field).onsets1, fTrialBehav1.(field).onsets2);
    excludedOffset = setdiff(fTrialBehav1.(field).offsets1, fTrialBehav1.(field).offsets2);
    for event = 1:numel(excludedOnset)
        onOnly = excludedOnset(event)-1;
        offOnly = excludedOffset(event)+1;
        excluded(count,:) = [offOnly-onOnly offOnly-onOnly+1];
        count = count +1;
    end
    IEItemp = fTrialBehav1.(field).onsets2(2:end)-fTrialBehav1.(field).offsets2(1:end-1); % for a better representation, best to use onsets/offsets2, instead of 1
    IEI = [IEI; IEItemp];
end

clear excludedOnset excludedOffset count t field on off IEItemp

statsBehSum = [];
statsBehSum(:,1)  = fR;
statsBehSum(:,2)  = sum(l);
statsBehSum(:,3)  = size(find(statsBeh(:,1)==0),1);
statsBehSum(:,4)  = size(find(statsBeh(:,1)==1),1);
statsBehSum(:,5)  = mean(statsBeh((statsBeh(:,1)==0),2));
statsBehSum(:,6)  = mean(statsBeh((statsBeh(:,1)==1),2));
statsBehSum(:,7)  = sum(statsBeh((statsBeh(:,1)==0),3));
statsBehSum(:,8)  = sum(statsBeh((statsBeh(:,1)==1),3));
statsBehSum(:,9)  = mean(statsBeh((statsBeh(:,1)==1),4));
statsBehSum(:,10) = max(statsBeh((statsBeh(:,1)==1),5));
statsBehSum(:,11) = eventCounter0-1;
statsBehSum(:,12) = eventCounter1-1;
statsBehSum(:,13) = size(excluded,1);
statsBehSum(:,14) = mean(excluded(:,1));
statsBehSum(:,15) = sum(excluded(:,2));
statsBehSum(:,16) = mean(IEI);

clear l excluded eventCounter0 eventCounter1

%% --- Save vars --- UGLYYY (consider eval?)
tbSpontAll = [];
tbSpontLedCl = [];
tbSpontLed = [];
if strcmp(groupName, 'spontAll') == 1
    tbSpontAll.group           = group;
    tbSpontAll.paramTrig       = paramTrig;
    tbSpontAll.movBinary       = movBinary;
    tbSpontAll.rates           = rates;
    tbSpontAll.RWhiskSpeed     = RWhiskSpeed;
    tbSpontAll.statsCa         = statsCa;
    tbSpontAll.statsCaSum      = statsCaSum;
    tbSpontAll.W               = W;
    tbSpontAll.WL              = WL;
    tbSpontAll.blockWS         = blockWS;
    tbSpontAll.blockWE         = blockWE;
    tbSpontAll.blockWLS        = blockWLS;
    tbSpontAll.blockWLE        = blockWLE;
    tbSpontAll.blockWLSresized = blockWLSresized;
    tbSpontAll.caAvW           = caAvW;
    tbSpontAll.caAvWL          = caAvWL;
    tbSpontAll.statsBeh        = statsBeh;
    tbSpontAll.statsBehSum     = statsBehSum;
    tbSpontAll.statsBehIEI     = IEI; %%
    tbSpontAll.corrWhisk1      = corrWhisk1;
    tbSpontAll.corrWhisk2      = corrWhisk2;
    tbSpontAll.rLR             = rLR;
    save([currExpVars.folder '\' baseSave 'trig' num2str(activity)], 'tbSpontAll')
end
if strcmp(groupName, 'spontLed') == 1
    tbSpontLed.group           = group;
    tbSpontLed.paramTrig       = paramTrig;
    tbSpontLed.movBinary       = movBinary;
    tbSpontLed.rates           = rates;
    tbSpontLed.RWhiskSpeed     = RWhiskSpeed;
    tbSpontLed.statsCa         = statsCa;
    tbSpontLed.statsCaSum      = statsCaSum;
    tbSpontLed.caAvW           = caAvW;
    tbSpontLed.caAvWL          = caAvWL;
    tbSpontLed.blockWS         = blockWS;
    tbSpontLed.blockWE         = blockWE;
    tbSpontLed.blockWLS        = blockWLS;
    tbSpontLed.blockWLE        = blockWLE;
    tbSpontLed.blockWLSresized = blockWLSresized;
    tbSpontLed.caAvW           = caAvW;
    tbSpontLed.caAvWL          = caAvWL;
    tbSpontLed.statsBeh        = statsBeh;
    tbSpontLed.statsBehSum     = statsBehSum;
    tbSpontLed.statsBehIEI     = IEI; %%
    tbSpontLed.corrWhisk1      = corrWhisk1;
    tbSpontLed.corrWhisk2      = corrWhisk2;
    tbSpontLed.rLR             = rLR;
    save([currExpVars.folder '\' baseSave 'trig' num2str(activity)], 'tbSpontLed', '-append')
end
if strcmp(groupName, 'spontLedCl') == 1
    tbSpontLedCl.group           = group;
    tbSpontLedCl.paramTrig       = paramTrig;
    tbSpontLedCl.movBinary       = movBinary;
    tbSpontLedCl.rates           = rates;
    tbSpontLedCl.RWhiskSpeed     = RWhiskSpeed;
    tbSpontLedCl.statsCa         = statsCa;
    tbSpontLedCl.statsCaSum      = statsCaSum;
    tbSpontLedCl.W               = W;
    tbSpontLedCl.WL              = WL;
    tbSpontLedCl.blockWS         = blockWS;
    tbSpontLedCl.blockWE         = blockWE;
    tbSpontLedCl.blockWLS        = blockWLS;
    tbSpontLedCl.blockWLE        = blockWLE;
    tbSpontLedCl.blockWLSresized = blockWLSresized;
    tbSpontLedCl.caAvW           = caAvW;
    tbSpontLedCl.caAvWL          = caAvWL;
    tbSpontLedCl.statsBeh        = statsBeh;
    tbSpontLedCl.statsBehSum     = statsBehSum;
    tbSpontLedCl.statsBehIEI     = IEI; %%
    tbSpontLedCl.corrWhisk1      = corrWhisk1;
    tbSpontLedCl.corrWhisk2      = corrWhisk2;
    tbSpontLedCl.rLR             = rLR;
    save([currExpVars.folder '\' baseSave 'trig' num2str(activity)], 'tbSpontLedCl', '-append')
end
% save([currExpVars.folder '\' baseSave 'trigBEH_' groupName], 'group', 'paramTrig', 'spksRate', 'RWhiskSpeed', 'statsCa', 'statsCaSum', 'blockWS', 'blockWE', 'blockWLS', 'blockWLE', 'blockWLSresized', 'statsBeh', 'statsBehSum', 'corrWhisk1', 'corrWhisk2', 'rLR')

%% --- Plotting ---
mkdir([currExpVars.folder '\Figs0\'])
mkdir([currExpVars.folder '\Figs1\'])

if ~isempty(blockWS)
    h1 = trigFig(blockWS, paramTrig, WCtrl(:,1), WDown(:,1), WUp(:,1), 0);
    saveas(h1, [currExpVars.folder '\Figs' num2str(activity) '\' baseSave '_tb_' groupName '_WOn'])
    
    h2 = trigFig(blockWE, paramTrig, WCtrl(:,1), WDown(:,1), WUp(:,1), 1);
    saveas(h2, [currExpVars.folder '\Figs' num2str(activity) '\' baseSave '_tb_' groupName '_WOff' ])
end

h3 = trigFig(blockWLS, paramTrig, WLCtrl(:,1), WLDown(:,1), WLUp(:,1), 0);
saveas(h3, [currExpVars.folder '\Figs' num2str(activity) '\' baseSave '_tb_' groupName '_WLOn'])

h4 = trigFig(blockWLE, paramTrig, WLCtrl(:,1), WLDown(:,1), WLUp(:,1), 2);
saveas(h4, [currExpVars.folder '\Figs' num2str(activity) '\' baseSave '_tb_' groupName '_WLOff' ])

h5 = trigFig(blockWLSresized, paramTrig, WLCtrl(:,1), WLDown(:,1), WLUp(:,1), 1);
saveas(h5, [currExpVars.folder '\Figs' num2str(activity) '\' baseSave '_tb_' groupName '_WLOnResized'])

close all
clear h1 h2 h3 h4 h5

%% Extract starter stats and figs
statsCaStarter = [];
if ~isempty(starterID)
    
    cell = find(goodCells==starterID);
    statsCaStarter = [statsCa(cell,:) rates(cell,:)];
    statsCaStarterSpont = statsCaStarter; % Really ugly; work it out
    save([currExpVars.folder '\' baseSave 'trig' num2str(activity)], 'statsCaStarterSpont', '-append')
    
    % Produce main triggered figures for the starter cell
    % Spont trig data
    if ~isempty(blockWS)
        h7 = trigFigCellByCell(blockWS, paramTrig, 0, cell);
        saveas(h7, [currExpVars.folder '\Figs' num2str(activity) '\' baseSave '_tb_' groupName '_Starter_WOn'])
        
        h8 = trigFigCellByCell(blockWE, paramTrig, 1, cell);
        saveas(h8, [currExpVars.folder '\Figs' num2str(activity) '\' baseSave '_tb_' groupName '_Starter_WOff'])
    end
    
    h9 = trigFigCellByCell(blockWLS, paramTrig, 0, cell);
    saveas(h9, [currExpVars.folder '\Figs' num2str(activity) '\' baseSave '_tb_' groupName '_Starter_WLOn'])
    
    h10 = trigFigCellByCell(blockWLE, paramTrig, 2, cell);
    saveas(h10, [currExpVars.folder '\Figs' num2str(activity) '\' baseSave '_tb_' groupName '_Starter_WLOff'])
    
    h11 = trigFigCellByCell(blockWLSresized, paramTrig, 1, cell);
    saveas(h11, [currExpVars.folder '\Figs' num2str(activity) '\' baseSave '_tb_' groupName '_Starter_WLOnResized'])
end
close all
clear statsCa h7 h8 h9 h10 h11
end
