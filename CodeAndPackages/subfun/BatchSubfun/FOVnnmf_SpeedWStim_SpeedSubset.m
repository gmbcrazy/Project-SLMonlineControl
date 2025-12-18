function FOVres = FOVnnmf_SpeedWStim_SpeedSubset(iFOV, SLMGroupTableTrial, NeuroTrace, BehTrace, iData, lowDimN, SponInd, ExcludeStimInd, SpeedPercentile)
    % FOVNNMF_SPEEDWSTIM_SPEEDSUBSET - Perform NNMF using high/low speed periods
    % 
    % Inputs:
    %   iFOV - Field of view ID
    %   SLMGroupTableTrial - Table with trial information
    %   NeuroTrace - Neural activity traces
    %   BehTrace - Behavioral traces (BehTrace{1} contains speed)
    %   iData - Data index
    %   lowDimN - Number of NNMF components
    %   SponInd - Spontaneous period indices
    %   ExcludeStimInd - Indices excluding stimulation
    %   SpeedPercentile - Percentage for high/low speed selection (e.g., 25 for top/bottom 25%)

    % Extract session-specific data
    tblFOV = SLMGroupTableTrial(SLMGroupTableTrial.Session == iFOV, :);
    TrialID = unique(tblFOV.TrialID);
    tbltemp = SLMGroupTableTrial(SLMGroupTableTrial.Session == iFOV & SLMGroupTableTrial.TrialID == TrialID(1), :);

    % Get speed data during spontaneous period
    SpeedSpon = mean(BehTrace(iFOV).Speed(SponInd,:),2);
    
    % Calculate speed percentiles
    lowSpeedThreshold = prctile(SpeedSpon(ExcludeStimInd), SpeedPercentile);
    highSpeedThreshold = prctile(SpeedSpon(ExcludeStimInd), 100 - SpeedPercentile);
    
    % Find low and high speed indices within SponInd
    lowSpeedIdx_rel = find(SpeedSpon <= lowSpeedThreshold);  % relative to SponInd
    highSpeedIdx_rel = find(SpeedSpon >= highSpeedThreshold);  % relative to SponInd
    
    % Convert to absolute indices
    lowSpeedIdx = SponInd(lowSpeedIdx_rel);
    highSpeedIdx = SponInd(highSpeedIdx_rel);
    
    % Combine low speed, high speed, and StimInd (union)
    % ExcludeStimInd should be intersected with valid data range
    trainingIdx = union(union(lowSpeedIdx, highSpeedIdx), setdiff(SponInd,ExcludeStimInd));
    trainingIdx = sort(trainingIdx);  % Sort for consistency
    
    % Prepare full data for reconstruction
    fullData = AmpNormalizeRow(NeuroTrace{iFOV}{iData}, [0 100]);
    
    % Prepare training data (subset based on speed)
    trainingData = fullData(:, trainingIdx);

    % NNMF for all neurons using training data
    [W, H_training] = nnmf(trainingData, lowDimN, 'algorithm', 'als', 'replicates', 20);

    % Reconstruct H for all time points using learned W
    H = W \ fullData;  % Solve W*H = fullData for H
    H(H < 0) = 0;  % Ensure non-negativity

    % Identify target and non-target cells
    I1 = find(tblFOV.TrialID == TrialID(1));
    NonTarget = find(tblFOV.NonTargetCell(I1) == 1);
    Target = find(tblFOV.TargetCell(I1) == 1);

  

    % Store results
    FOVres.W = W;
    FOVres.H = H;
    FOVres.H_training = H_training;
    FOVres.trainingIdx = trainingIdx;
    FOVres.Nontarget = NonTarget;
    FOVres.Target = Target;
    FOVres.SpeedPercentile = SpeedPercentile;
    FOVres.lowSpeedIdx = lowSpeedIdx;
    FOVres.highSpeedIdx = highSpeedIdx;

    % NNMF for only non-target neurons using training data
    SubData_training = trainingData(NonTarget, :);
    [W_NonTarget, H_NonTarget_training] = nnmf(SubData_training, lowDimN, 'algorithm', 'als', 'replicates', 20);
    
    % Reconstruct H_NonTarget for all time points
    SubData_full = fullData(NonTarget, :);
    H_NonTarget = W_NonTarget \ SubData_full;
    H_NonTarget(H_NonTarget < 0) = 0;
    
    FOVres.W_NonTarget = W_NonTarget;
    FOVres.H_NonTarget = H_NonTarget;
    FOVres.H_NonTarget_training = H_NonTarget_training;

    % Compute correlation with behavior (speed and stim)
    HH = H';
    SpeedSpon_full = nanmean(BehTrace(iFOV).Speed(ExcludeStimInd, :), 2);
    rSpeed = corr(HH(ExcludeStimInd, 1:lowDimN), SpeedSpon_full, 'rows', 'complete');

    maxLag = 10;
    for PCI = 1:lowDimN
        [c(:, PCI), lags] = xcorr(AmpNormalize(HH(SponInd, PCI), [0 100]), ...
            AmpNormalize(nanmean(BehTrace(iFOV).Stim(SponInd, :), 2), [0 100]), maxLag, 'coeff');
        PostI = find(lags >= 0 & lags <= 10);
        PreI = find(lags < 0);
        [~, i1] = max(abs(c(PostI, PCI)) - mean(c(PreI, PCI)));
        rStim(PCI, 1) = c(PostI(i1), PCI) - mean(AmpNormalize(HH(SponInd, PCI)));
    end

    [~, pcaIspeed] = max(rSpeed);

    [~, pcaIstim] = sort(rStim, 'descend');
    pcaIstim(pcaIstim == pcaIspeed) = [];
    pcaIstim = pcaIstim(1);

    FOVres.SpeedScore = rSpeed([pcaIspeed pcaIstim]);
    FOVres.SensoryScore = rStim([pcaIspeed pcaIstim]);

    % Repeat for non-target only (NonTarget data)
    HH_NonTarget = H_NonTarget';
    rSpeed_NonTarget = corr(HH_NonTarget(ExcludeStimInd, 1:lowDimN), SpeedSpon_full, 'rows', 'complete');

    for PCI = 1:lowDimN
        [c_NonTarget(:, PCI), lags] = xcorr(AmpNormalize(HH_NonTarget(SponInd, PCI), [0 100]), ...
            AmpNormalize(nanmean(BehTrace(iFOV).Stim(SponInd, :), 2), [0 100]), maxLag, 'coeff');
        PostI = find(lags >= 0 & lags <= 10);
        PreI = find(lags < 0);
        [~, i1] = max(abs(c_NonTarget(PostI, PCI)) - mean(c_NonTarget(PreI, PCI)));
        rStim_NonTarget(PCI, 1) = c_NonTarget(PostI(i1), PCI) - mean(AmpNormalize(HH_NonTarget(SponInd, PCI)));
    end

    [~, pcaIspeed_NonTarget] = max(rSpeed_NonTarget);
    [~, pcaIstim_NonTarget] = sort(rStim_NonTarget, 'descend');
    pcaIstim_NonTarget(pcaIstim_NonTarget == pcaIspeed_NonTarget) = [];
    pcaIstim_NonTarget = pcaIstim_NonTarget(1);

    FOVres.SpeedScore_NonTarget = rSpeed_NonTarget([pcaIspeed_NonTarget pcaIstim_NonTarget]);
    FOVres.SensoryScore_NonTarget = rStim_NonTarget([pcaIspeed_NonTarget pcaIstim_NonTarget]);

    % Linear regression with speed
    [FOVres.SpeedReg.B, FOVres.SpeedReg.BINT, FOVres.SpeedReg.R, FOVres.SpeedReg.RINT, FOVres.SpeedReg.STATS] = ...
        regress(H(pcaIspeed, ExcludeStimInd)', [ones(length(SpeedSpon_full), 1) SpeedSpon_full]);

    % [FOVres.SpeedRegNontarget.B, FOVres.SpeedRegNontarget.BINT, FOVres.SpeedRegNontarget.R, FOVres.SpeedRegNontarget.RINT, FOVres.SpeedRegNontarget.STATS] = ...
    %     regress(HNontarget(pcaIspeed, ExcludeStimInd)', [ones(length(SpeedSpon_full), 1) SpeedSpon_full]);

    [FOVres.SpeedReg_NonTarget.B, FOVres.SpeedReg_NonTarget.BINT, FOVres.SpeedReg_NonTarget.R, FOVres.SpeedReg_NonTarget.RINT, FOVres.SpeedReg_NonTarget.STATS] = ...
        regress(H_NonTarget(pcaIspeed_NonTarget, ExcludeStimInd)', [ones(length(SpeedSpon_full), 1) SpeedSpon_full]);

    % Store indices
    FOVres.ComI = [pcaIspeed pcaIstim];
    FOVres.ComI_NonTarget = [pcaIspeed_NonTarget pcaIstim_NonTarget];

    FOVres.FOVid = iFOV;

end