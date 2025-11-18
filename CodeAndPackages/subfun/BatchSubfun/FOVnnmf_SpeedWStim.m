function FOVres = FOVnnmf_SpeedWStim(iFOV, SLMGroupTableTrial, NeuroTrace, BehTrace, iData, lowDimN, SponInd, ExcludeStimInd)
    % PROCESSFOVITERATION - Perform dimension reduction and compute correlations
    % between neural and behavioral data for a given field of view (FOV).

    % Extract session-specific data
    tblFOV = SLMGroupTableTrial(SLMGroupTableTrial.Session == iFOV, :);
    TrialID = unique(tblFOV.TrialID);
    tbltemp = SLMGroupTableTrial(SLMGroupTableTrial.Session == iFOV & SLMGroupTableTrial.TrialID == TrialID(1), :);

    % NNMF for all neurons
    traningData=AmpNormalizeRow(NeuroTrace{iFOV}{iData}(:, SponInd),[0 100]);


    [W, H] = nnmf(traningData, lowDimN, 'algorithm', 'als', 'replicates', 20);

    % Identify target and non-target cells
    I1 = find(tblFOV.TrialID == TrialID(1));
    NonTarget = find(tblFOV.NonTargetCell(I1) == 1);
    Target = find(tblFOV.TargetCell(I1) == 1);

    % Compute NNMF components for non-target cells
    WNontarget = W(NonTarget, :);
    pinvW = pinv(W);
    subNonTargetData = traningData;
    subNonTargetData(Target, :) = 0;
    HNontarget = (pinvW * subNonTargetData);

    % Store results
    FOVres.W = W;
    FOVres.H = H;
    FOVres.WNontarget = WNontarget;
    FOVres.HNontarget = HNontarget;
    FOVres.Nontarget = NonTarget;
    FOVres.Target = Target;

    % NNMF for only non-target neurons
    SubData = traningData(NonTarget, :);


    [W_NonTarget, H_NonTarget] = nnmf(SubData, lowDimN, 'algorithm', 'als', 'replicates', 20);
    FOVres.W_NonTarget = W_NonTarget;
    FOVres.H_NonTarget = H_NonTarget;

    % Compute correlation with behavior (speed and stim)
    HH = H';
    SpeedSpon = nanmean(BehTrace(iFOV).Speed(ExcludeStimInd, :), 2);
    rSpeed = corr(HH(ExcludeStimInd, 1:lowDimN), SpeedSpon, 'rows', 'complete');

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
    [~, pcaIstim] = max(rStim);

    FOVres.SpeedScore = rSpeed([pcaIspeed pcaIstim]);
    FOVres.SensoryScore = rStim([pcaIspeed pcaIstim]);

    % Repeat for non-target only (NonTarget data)
    HH_NonTarget = H_NonTarget';
    rSpeed_NonTarget = corr(HH_NonTarget(ExcludeStimInd, 1:lowDimN), SpeedSpon, 'rows', 'complete');

    for PCI = 1:lowDimN
        [c_NonTarget(:, PCI), lags] = xcorr(AmpNormalize(HH_NonTarget(SponInd, PCI), [0 100]), ...
            AmpNormalize(nanmean(BehTrace(iFOV).Stim(SponInd, :), 2), [0 100]), maxLag, 'coeff');
        PostI = find(lags >= 0 & lags <= 10);
        PreI = find(lags < 0);
        [~, i1] = max(abs(c_NonTarget(PostI, PCI)) - mean(c_NonTarget(PreI, PCI)));
        rStim_NonTarget(PCI, 1) = c_NonTarget(PostI(i1), PCI) - mean(AmpNormalize(HH_NonTarget(SponInd, PCI)));
    end

    [~, pcaIspeed_NonTarget] = max(rSpeed_NonTarget);
    [~, pcaIstim_NonTarget] = max(rStim_NonTarget);

    FOVres.SpeedScore_NonTarget = rSpeed_NonTarget([pcaIspeed_NonTarget pcaIstim_NonTarget]);
    FOVres.SensoryScore_NonTarget = rStim_NonTarget([pcaIspeed_NonTarget pcaIstim_NonTarget]);

    % Linear regression with speed
    [FOVres.SpeedReg.B, FOVres.SpeedReg.BINT, FOVres.SpeedReg.R, FOVres.SpeedReg.RINT, FOVres.SpeedReg.STATS] = ...
        regress(H(pcaIspeed, ExcludeStimInd)', [ones(length(SpeedSpon), 1) SpeedSpon]);

    [FOVres.SpeedRegNontarget.B, FOVres.SpeedRegNontarget.BINT, FOVres.SpeedRegNontarget.R, FOVres.SpeedRegNontarget.RINT, FOVres.SpeedRegNontarget.STATS] = ...
        regress(HNontarget(pcaIspeed, ExcludeStimInd)', [ones(length(SpeedSpon), 1) SpeedSpon]);

    [FOVres.SpeedReg_NonTarget.B, FOVres.SpeedReg_NonTarget.BINT, FOVres.SpeedReg_NonTarget.R, FOVres.SpeedReg_NonTarget.RINT, FOVres.SpeedReg_NonTarget.STATS] = ...
        regress(H_NonTarget(pcaIspeed_NonTarget, ExcludeStimInd)', [ones(length(SpeedSpon), 1) SpeedSpon]);

    % Store indices
    FOVres.ComI = [pcaIspeed pcaIstim];
    FOVres.ComI_NonTarget = [pcaIspeed_NonTarget pcaIstim_NonTarget];


    FOVres.FOVid=iFOV;

end


% function plotFOVResults(FOVres, iFOV, SponInd, BehTrace, ProcessPar, SaveFunSub, NDataName, iData, IndexFOVNeed)
%     % PLOTFOVRESULTS - Visualization separated from data processing logic
%     % (implement plotting from original loop here)
% end
