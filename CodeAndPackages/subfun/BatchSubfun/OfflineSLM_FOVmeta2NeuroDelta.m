function [statGroupRes, CellIDFOV] = OfflineSLM_FOVmeta2NeuroDelta(InputFOVData, Par, PSTHparam)
% STRATEGY A
% Compute per-cell Δ (post - pre) for each stimulation group within each
% (VolOut, AwakeState) condition, plus a zero-power control column, and a
% normalized delta (group minus zero-power).
%
% Required fields (in InputFOVData / Par):
%   InputFOVData.AlignedNData{fov}         : [cells x time x trials]
%   InputFOVData.AlignedInfoTableFOV{fov}  : table with columns Group, PowerZero, VolOut, AwakeState, PointTargetCell, ...
%   InputFOVData.AlignedInfoTable          : global table with column iFOV (trial -> FOV id)
%   InputFOVData.AlignedSpeedMeta          : [time x allTrials]
%   InputFOVData.NeuroPos3DMeta            : [:,4] contains FOV id
%   Par.GroupList, Par.VolOut, Par.AwakeState, Par.PostPreDiffSpeedTh
%   PSTHparam.PreSLMCal, .PreTestFrame, .TestStepFrame
%
% Outputs:
%   statGroupRes(iVol,iAwake) with fields:
%       p, delta, deltaNorm, GroupResponse, GroupSpeed, GroupSampleN,
%       GroupSpeedTrial (concatenated [preSpeed postSpeed] per-trial rows),
%       speeddelta (per-FOV speed delta, replicated to its cells)
%   CellIDFOV : [NCell x 2] = [withinFOVindex, FOVid]

StructToVars(Par);
StructToVars(InputFOVData);

nGroupsReal   = numel(GroupList);
FunNum        = nGroupsReal + 1;      % +1 for zero-power
NCell         = size(NeuroPos3DMeta,1);

% Time windows (exactly the same indices Strategy-B Win1 should use)
preIdx  = (PSTHparam.PreSLMCal - PSTHparam.PreTestFrame + 1) : PSTHparam.PreSLMCal;
postIdx = PSTHparam.PreSLMCal + (1:PSTHparam.TestStepFrame);

% FOV bookkeeping
nFOV = numel(AlignedNData);
FOVCellIdx  = cell(nFOV,1);
FOVEventIdx = cell(nFOV,1);
for iFOV = 1:nFOV
    FOVCellIdx{iFOV}  = find(NeuroPos3DMeta(:,4) == iFOV);
    FOVEventIdx{iFOV} = find(AlignedInfoTable.iFOV == iFOV);
end

% Map each neuron to its within-FOV index and FOV id (once)
CellIDFOV = zeros(NCell,2);
for iFOV = 1:nFOV
    iFOVneuro = FOVCellIdx{iFOV};
    CellIDFOV(iFOVneuro,1) = (1:numel(iFOVneuro)).';
    CellIDFOV(iFOVneuro,2) = iFOV;
end

% Length of PSTH axis (for storage)
nTime = PSTHparam.PreSLMCal + PSTHparam.PostSLMCal;

for iVol = 1:numel(VolOut)
    for iAw = 1:numel(AwakeState)

        % Allocate outputs for this condition
        statGroupRes(iVol,iAw).p             = ones(NCell, FunNum);
        statGroupRes(iVol,iAw).delta         = zeros(NCell, FunNum);
        statGroupRes(iVol,iAw).deltaNorm     = zeros(NCell, nGroupsReal);
        statGroupRes(iVol,iAw).GroupSampleN  = zeros(NCell, FunNum);
        statGroupRes(iVol,iAw).GroupResponse = nan(NCell, nTime, FunNum);
        statGroupRes(iVol,iAw).GroupSpeed    = nan(NCell, nTime, FunNum);
        statGroupRes(iVol,iAw).GroupSpeedTrial = cell(FunNum,1);
        statGroupRes(iVol,iAw).speeddelta    = zeros(NCell, FunNum);

        % ---------- Real groups (PowerZero==0) ----------
        for iFun = 1:nGroupsReal
            SpeedGroupAcc = [];

            for iFOV = 1:nFOV
                T = AlignedInfoTableFOV{iFOV};
                I1 = find(T.Group==GroupList(iFun) & T.PowerZero==0 & ...
                          T.VolOut==VolOut(iVol)   & T.AwakeState==AwakeState(iAw));
                if isempty(I1), continue; end

                iFOVneuro = FOVCellIdx{iFOV};
                iFOVEvent = FOVEventIdx{iFOV};

                % FOV-scoped speed
                SpeedFOV = AlignedSpeedMeta(:, iFOVEvent);
                preSp  = mean(SpeedFOV(preIdx, I1), 1, 'omitnan');
                postSp = mean(SpeedFOV(postIdx, I1), 1, 'omitnan');

                % Speed gating only when awake
                if AwakeState(iAw) == 1
                    keep = (postSp - preSp) <= PostPreDiffSpeedTh;
                    I1    = I1(keep);
                    preSp = preSp(keep);
                    postSp= postSp(keep);
                end
                if isempty(I1), continue; end

                % Accumulate [pre post] per trial for this group
                SpeedGroupAcc = [SpeedGroupAcc; [preSp(:) postSp(:)]]; %#ok<AGROW>

                % FOV-level speed delta
                SpeedDelta = mean(postSp) - mean(preSp);
                statGroupRes(iVol,iAw).speeddelta(iFOVneuro, iFun) = SpeedDelta;

                % Neural responses
                X = AlignedNData{iFOV}; % [cells x time x trials]

                if numel(I1) == 1
                    % Single trial: replicate speed trace over cells; unpaired test across frames
                    statGroupRes(iVol,iAw).GroupSpeed(iFOVneuro,:,iFun)    = ...
                        repmat(SpeedFOV(:, I1).', numel(iFOVneuro), 1);
                    statGroupRes(iVol,iAw).GroupResponse(iFOVneuro,:,iFun) = ...
                        squeeze(X(:,:,I1));

                    statGroupRes(iVol,iAw).GroupSampleN(iFOVneuro,iFun)   = 1;

                    preDat  = squeeze(X(:, preIdx,  I1));   % [cells x preFrames]
                    postDat = squeeze(X(:, postIdx, I1));   % [cells x postFrames]

                    pVec     = nan(numel(iFOVneuro),1);
                    dVec     = nan(numel(iFOVneuro),1);
                    for j = 1:numel(iFOVneuro)
                        [~, pVec(j)] = ttest2(postDat(j,:).', preDat(j,:).'); % unpaired
                        dVec(j)      = mean(postDat(j,:), 'omitnan') - mean(preDat(j,:), 'omitnan');
                    end
                    statGroupRes(iVol,iAw).p(iFOVneuro,    iFun) = pVec;
                    statGroupRes(iVol,iAw).delta(iFOVneuro,iFun) = dVec;

                else
                    % Multi-trial: average within windows per trial, paired across trials
                    statGroupRes(iVol,iAw).GroupSpeed(iFOVneuro,:,iFun)    = ...
                        repmat(mean(SpeedFOV(:, I1), 2, 'omitnan').', numel(iFOVneuro), 1);
                    statGroupRes(iVol,iAw).GroupResponse(iFOVneuro,:,iFun) = ...
                        mean(X(:,:,I1), 3, 'omitnan');

                    statGroupRes(iVol,iAw).GroupSampleN(iFOVneuro,iFun)   = numel(I1);

                    preAvg  = squeeze(mean(X(:, preIdx,  I1), 2, 'omitnan')); % cells x trials
                    postAvg = squeeze(mean(X(:, postIdx, I1), 2, 'omitnan'));
                    if isvector(preAvg), preAvg = preAvg(:); postAvg = postAvg(:); end

                    pVec = nan(numel(iFOVneuro),1);
                    dVec = nan(numel(iFOVneuro),1);
                    for j = 1:numel(iFOVneuro)
                        [~, pVec(j)] = ttest(postAvg(j,:), preAvg(j,:));      % paired across trials
                        dVec(j)      = mean(postAvg(j,:) - preAvg(j,:), 'omitnan');
                    end
                    statGroupRes(iVol,iAw).p(iFOVneuro,    iFun) = pVec;
                    statGroupRes(iVol,iAw).delta(iFOVneuro,iFun) = dVec;
                end
            end

            statGroupRes(iVol,iAw).GroupSpeedTrial{iFun} = SpeedGroupAcc;
        end

        % ---------- Zero-power control as last column ----------
        iFunZ = nGroupsReal + 1;
        SpeedGroupAcc = [];

        for iFOV = 1:nFOV
            T = AlignedInfoTableFOV{iFOV};
            I1 = find(T.PowerZero==1 & T.VolOut==VolOut(iVol) & T.AwakeState==AwakeState(iAw));
            if isempty(I1), continue; end

            iFOVneuro = FOVCellIdx{iFOV};
            iFOVEvent = FOVEventIdx{iFOV};

            SpeedFOV = AlignedSpeedMeta(:, iFOVEvent);
            preSp  = mean(SpeedFOV(preIdx, I1), 1, 'omitnan');
            postSp = mean(SpeedFOV(postIdx, I1), 1, 'omitnan');

            % Speed gating only when awake (keep all sham otherwise if anesthetized)
            if AwakeState(iAw) == 1
                keep = (postSp - preSp) <= PostPreDiffSpeedTh;
                I1    = I1(keep);
                preSp = preSp(keep);
                postSp= postSp(keep);
            end
            if isempty(I1), continue; end

            SpeedGroupAcc = [SpeedGroupAcc; [preSp(:) postSp(:)]]; %#ok<AGROW>

            SpeedDelta = mean(postSp) - mean(preSp);
            statGroupRes(iVol,iAw).speeddelta(iFOVneuro, iFunZ) = SpeedDelta;

            X = AlignedNData{iFOV};

            if numel(I1) == 1
                statGroupRes(iVol,iAw).GroupSpeed(iFOVneuro,:,iFunZ)    = ...
                    repmat(SpeedFOV(:, I1).', numel(iFOVneuro), 1);
                statGroupRes(iVol,iAw).GroupResponse(iFOVneuro,:,iFunZ) = ...
                    squeeze(X(:,:,I1));

                statGroupRes(iVol,iAw).GroupSampleN(iFOVneuro,iFunZ)   = 1;

                preDat  = squeeze(X(:, preIdx,  I1));
                postDat = squeeze(X(:, postIdx, I1));

                pVec = nan(numel(iFOVneuro),1);
                dVec = nan(numel(iFOVneuro),1);
                for j = 1:numel(iFOVneuro)
                    [~, pVec(j)] = ttest2(postDat(j,:).', preDat(j,:).');
                    dVec(j)      = mean(postDat(j,:), 'omitnan') - mean(preDat(j,:), 'omitnan');
                end
                statGroupRes(iVol,iAw).p(iFOVneuro,    iFunZ) = pVec;
                statGroupRes(iVol,iAw).delta(iFOVneuro,iFunZ) = dVec;

            else
                statGroupRes(iVol,iAw).GroupResponse(iFOVneuro,:,iFunZ) = ...
                    mean(X(:,:,I1), 3, 'omitnan');
                statGroupRes(iVol,iAw).GroupSpeed(iFOVneuro,:,iFunZ)    = ...
                    repmat(mean(SpeedFOV(:, I1), 2, 'omitnan').', numel(iFOVneuro), 1);
                statGroupRes(iVol,iAw).GroupSampleN(iFOVneuro,iFunZ)    = numel(I1);

                preAvg  = squeeze(mean(X(:, preIdx,  I1), 2, 'omitnan'));
                postAvg = squeeze(mean(X(:, postIdx, I1), 2, 'omitnan'));
                if isvector(preAvg), preAvg = preAvg(:); postAvg = postAvg(:); end

                pVec = nan(numel(iFOVneuro),1);
                dVec = nan(numel(iFOVneuro),1);
                for j = 1:numel(iFOVneuro)
                    [~, pVec(j)] = ttest(postAvg(j,:), preAvg(j,:));
                    dVec(j)      = mean(postAvg(j,:) - preAvg(j,:), 'omitnan');
                end
                statGroupRes(iVol,iAw).p(iFOVneuro,    iFunZ) = pVec;
                statGroupRes(iVol,iAw).delta(iFOVneuro,iFunZ) = dVec;
            end
        end

        statGroupRes(iVol,iAw).GroupSpeedTrial{iFunZ} = SpeedGroupAcc;

        % Normalized deltas: each real group minus the zero-power column
        statGroupRes(iVol,iAw).deltaNorm = ...
            statGroupRes(iVol,iAw).delta(:,1:nGroupsReal) - ...
            statGroupRes(iVol,iAw).delta(:,end);
    end
end
end


%% Originalversion, Lu Zhang, Oct 21 2025
% function [statGroupRes,CellIDFOV]=OfflineSLM_FOVmeta2NeuroDelta(InputFOVData, Par,PSTHparam)
% 
% StructToVars(Par);
% StructToVars(InputFOVData);
% 
% GroupMetaName = [GroupLabel {'FakeSLM'}];
% GroupMetaColor = [GroupColor; PowerZeroColor];
% NCell = size(InputFOVData.NeuroPos3DMeta,1);
% FunNum = length(GroupList) + 1;
% TimBinFrame = -PSTHparam.PreSLMCal:PSTHparam.PostSLMCal-1;
% 
% TestStepFrame=PSTHparam.TestStepFrame;
%         % clear statGroupRes
% CellIDFOV=zeros(NCell,2);
% 
% for iVol = 1:length(VolOut)
%     for iAwake = 1:length(AwakeState)
% 
%         % clear GroupResponse GroupSpeed statGroupRes GroupSampleN
% 
%         statGroupRes(iVol,iAwake).p         = ones(NCell, FunNum);
%         statGroupRes(iVol,iAwake).delta     = zeros(NCell, FunNum);
%         statGroupRes(iVol,iAwake).speeddelta= zeros(NCell, FunNum);
%         statGroupRes(iVol,iAwake).GroupSampleN = zeros(NCell, FunNum);
%         statGroupRes(iVol,iAwake).GroupSpeed = nan(NCell, length(TimBinFrame), FunNum);
%         statGroupRes(iVol,iAwake).GroupResponse = nan(NCell, length(TimBinFrame), FunNum);
%         statGroupRes(iVol,iAwake).GroupSpeedTrial = {};
% 
% 
%        for iFun = 1:length(GroupList)
%            I1All=find(AlignedInfoTable.Group==GroupList(iFun)&AlignedInfoTable.PowerZero==0&AlignedInfoTable.VolOut==VolOut(iVol)&AlignedInfoTable.AwakeState==AwakeState(iAwake));
%             % GroupSampleN(iFun)=length(I1);
% 
%            preSLMSpeed=mean(squeeze(AlignedSpeedMeta((PSTHparam.PreSLMCal-PSTHparam.PreTestFrame+1):PSTHparam.PreSLMCal,I1All)),1);
%            postSLMSpeed=mean(squeeze(AlignedSpeedMeta(PSTHparam.PreSLMCal+[1:TestStepFrame],I1All)),1);
% 
%            % % if iAwake==1
%            % % I1All=I1All(find(abs(postSLMSpeed-preSLMSpeed)<=PostPreDiffSpeedTh));
%            % % end
%            % % 
%            % % GroupSpeed{iFun}=AlignedSpeedMeta(:,I1All)';
%            % % GroupSampleN(iFun)=length(I1All);
%            % iFun
%            SpeedGroup{iFun}=[];
%            for iFOV = 1:length(InputFOVData.AlignedNData)
%                % iFOV
%                 AlignedInfoTableFOV=InputFOVData.AlignedInfoTableFOV{iFOV};
%                 I1=find(AlignedInfoTableFOV.Group==GroupList(iFun)&AlignedInfoTableFOV.PowerZero==0&AlignedInfoTableFOV.VolOut==VolOut(iVol)&AlignedInfoTableFOV.AwakeState==AwakeState(iAwake));
%                 AlignedtempNData = InputFOVData.AlignedNData{iFOV};
%                 iFOVneuro=find(InputFOVData.NeuroPos3DMeta(:,4)==iFOV);
%                 iFOVEvent = find(AlignedInfoTable.iFOV==iFOV);
% 
%                 if iFun==1
%                     CellIDFOV(iFOVneuro,1)=[1:length(iFOVneuro)]';
%                     CellIDFOV(iFOVneuro,2)=iFOV;
%                 end
%                 AlignedSpeedFOV=AlignedSpeedMeta(:,iFOVEvent);
%                 preSLMSpeed=mean(squeeze(AlignedSpeedFOV((PSTHparam.PreSLMCal-PSTHparam.PreTestFrame+1):PSTHparam.PreSLMCal,I1)),1);
%                 postSLMSpeed=mean(squeeze(AlignedSpeedFOV(PSTHparam.PreSLMCal+[1:TestStepFrame],I1)),1);
% 
% 
%                 if iAwake==1
%                    I1temp=find(abs(postSLMSpeed-preSLMSpeed)<=PostPreDiffSpeedTh);
%                    I1=I1(I1temp);
%                 end
%                 preSLMSpeed=preSLMSpeed(I1temp);
%                 postSLMSpeed=postSLMSpeed(I1temp);
%                 SpeedGroup{iFun}=[SpeedGroup{iFun};[preSLMSpeed(:) postSLMSpeed(:)]];
% 
% 
% 
% 
%                 [~,pSpeed(iFun, iFOV),~,statSpeed(iFun, iFOV)]=ttest(postSLMSpeed,preSLMSpeed);
%                 SpeedDelta = mean(postSLMSpeed) - mean(preSLMSpeed);
%                 % size(statGroupRes(iVol,iAwake).speeddelta)
%                 statGroupRes(iVol,iAwake).speeddelta(iFOVneuro,iFun)=SpeedDelta;
% 
%                 if length(I1)==1
%                    % statGroupRes(iVol,iAwake).GroupResponse(iFOVneuro,:,iFun)=squeeze(AlignedtempNData(:,:,I1));
%                    statGroupRes(iVol,iAwake).GroupSpeed(iFOVneuro,:,iFun)=repmat(AlignedSpeedMeta(:,I1)',length(iFOVneuro), 1);
%                    statGroupRes(iVol,iAwake).GroupResponse(iFOVneuro,:,iFun)=squeeze(AlignedtempNData(:,:,I1));
%                    statGroupRes(iVol,iAwake).GroupSampleN(iFOVneuro,iFun)=length(I1);
% 
%                    preSLMdata=squeeze(AlignedtempNData(:,(PSTHparam.PreSLMCal-PSTHparam.PreTestFrame+1):PSTHparam.PreSLMCal,I1));
%                    postSLMdata=squeeze(AlignedtempNData(:,PSTHparam.PreSLMCal+[1:TestStepFrame],I1));
%                    for jCell=1:sum(InputFOVData.NeuroPos3DMeta(:,4)==iFOV)
%                        temp1=preSLMdata(jCell,:,:);
%                        temp2=postSLMdata(jCell,:,:);
%                        [~,p(jCell,1),~,t(jCell)]=ttest(temp2(:),temp1(:));
%                        temp3(jCell)=mean(temp2(:))-mean(temp1(:));
%                    end
%                 elseif length(I1)>1
%                    % GroupSpeed{iFun}=AlignedSpeedMeta(:,I1)';
%                    statGroupRes(iVol,iAwake).GroupSpeed(iFOVneuro,:,iFun)=repmat(nanmean(AlignedSpeedMeta(:,I1),2)',length(iFOVneuro), 1);
%                    statGroupRes(iVol,iAwake).GroupResponse(iFOVneuro,:,iFun)=nanmean(AlignedtempNData(:,:,I1),3);  
%                    preSLMdata=squeeze(AlignedtempNData(:,(PSTHparam.PreSLMCal-PSTHparam.PreTestFrame+1):PSTHparam.PreSLMCal,I1));
%                    postSLMdata=squeeze(AlignedtempNData(:,PSTHparam.PreSLMCal+[1:TestStepFrame],I1));
%                    for jCell=1:sum(InputFOVData.NeuroPos3DMeta(:,4)==iFOV)
%                        temp1=squeeze(nanmean(preSLMdata(jCell,:,:),2));
%                        temp2=squeeze(nanmean(postSLMdata(jCell,:,:),2));
%                        [~,p(jCell,1),~,t(jCell)]=ttest(temp2(:),temp1(:));
%                        temp3(jCell)=mean(temp2(:))-mean(temp1(:));
%                    end
%                    statGroupRes(iVol,iAwake).p(iFOVneuro,iFun)=p;
%                    % statGroupRes(iFun).t=t;
%                    statGroupRes(iVol,iAwake).delta(iFOVneuro,iFun)=temp3;
%                    statGroupRes(iVol,iAwake).GroupSampleN(iFOVneuro,iFun)=length(I1);
%                    clear p t temp3;
% 
% 
%                 else
%                    % statGroupRes(iFun).p=zeros(size(iscell))+1;
%                    % statGroupRes(iFun).delta=zeros(size(iscell));
% 
%                 end
% 
%         end
% 
% 
% 
%         end
% 
% 
% 
%         % statGroupRes(iFun).p=zeros(NCell,1)+1;
%         % statGroupRes(iFun).delta=zeros(NCell,1);
%         % GroupResponse{iFun}=zeros(NCell,length(TimBinFrame));
% 
%         iFun=iFun+1;
% %%Add Zero power as addtional group
%         I1All=find(AlignedInfoTable.PowerZero==1&AlignedInfoTable.VolOut==VolOut(iVol)&AlignedInfoTable.AwakeState==AwakeState(iAwake));
%         preSLMSpeed=mean(squeeze(AlignedSpeedMeta((PSTHparam.PreSLMCal-PSTHparam.PreTestFrame+1):PSTHparam.PreSLMCal,I1All)),1);
%         postSLMSpeed=mean(squeeze(AlignedSpeedMeta(PSTHparam.PreSLMCal+[1:TestStepFrame],I1All)),1);
%         % if iAwake==1
%         %    I1All=I1All(find(abs(postSLMSpeed-preSLMSpeed)<=PostPreDiffSpeedTh));
%         % end
%         % GroupSpeed{iFun}=AlignedSpeedMeta(:,I1All)';
% 
%         % GroupSampleN(iFun)=length(I1All);
%         SpeedGroup{iFun}=[];
%         for iFOV = 1:length(InputFOVData.AlignedNData)
%                 AlignedInfoTableFOV=InputFOVData.AlignedInfoTableFOV{iFOV};
%                 I1=find(AlignedInfoTableFOV.PowerZero==1&AlignedInfoTableFOV.VolOut==VolOut(iVol)&AlignedInfoTableFOV.AwakeState==AwakeState(iAwake));
%                 AlignedtempNData = InputFOVData.AlignedNData{iFOV};
%                 iFOVneuro=find(InputFOVData.NeuroPos3DMeta(:,4)==iFOV);
% 
% 
%                 iFOVEvent = find(AlignedInfoTable.iFOV==iFOV);
%                 AlignedSpeedFOV=AlignedSpeedMeta(:,iFOVEvent);
%                 preSLMSpeed=mean(squeeze(AlignedSpeedFOV((PSTHparam.PreSLMCal-PSTHparam.PreTestFrame+1):PSTHparam.PreSLMCal,I1)),1);
%                 postSLMSpeed=mean(squeeze(AlignedSpeedFOV(PSTHparam.PreSLMCal+[1:TestStepFrame],I1)),1);
% 
% 
% 
%                 if iAwake==1
%                     I1temp=find(abs(postSLMSpeed-preSLMSpeed)<=PostPreDiffSpeedTh);
%                     I1=I1(I1temp);
%                 end
%                 preSLMSpeed=preSLMSpeed(I1temp);
%                 postSLMSpeed=postSLMSpeed(I1temp);
%                 SpeedGroup{iFun}=[SpeedGroup{iFun};[preSLMSpeed(:) postSLMSpeed(:)]];
% 
% 
% 
%                 [~,pSpeed(iFun,iFOV),~,statSpeed(iFun,iFOV)]=ttest(postSLMSpeed,preSLMSpeed);
%                 SpeedDelta = mean(postSLMSpeed) - mean(preSLMSpeed);
%                 statGroupRes(iVol,iAwake).speeddelta(iFOVneuro,iFun)=SpeedDelta;
% 
% 
% 
% 
% 
% 
%             if length(I1)==1
%                % GroupResponse{iFun}=squeeze(AlignedtempNData(:,:,I1));
%                % GroupSpeed{iFun}=AlignedSpeedMeta(:,I1)';
%                statGroupRes(iVol,iAwake).GroupSpeed(iFOVneuro,:,iFun)=repmat(AlignedSpeedMeta(:,I1)',length(iFOVneuro),1);
%                statGroupRes(iVol,iAwake).GroupResponse(iFOVneuro,:,iFun)=squeeze(AlignedtempNData(:,:,I1));
%                statGroupRes(iVol,iAwake).GroupSampleN(iFOVneuro,iFun)=length(I1);
% 
%             elseif length(I1)>1
%                % GroupResponse{iFun}=nanmean(AlignedtempNData(:,:,I1),3);    
%                statGroupRes(iVol,iAwake).GroupResponse(iFOVneuro,:,iFun)=nanmean(AlignedtempNData(:,:,I1),3);  
%                % GroupSpeed{iFun}=AlignedSpeedMeta(:,I1)';
%                % statGroupRes(iVol,iAwake).GroupSpeed(iFOVneuro,:,iFun)=nanmean(AlignedSpeedMeta(:,I1),2);
%                statGroupRes(iVol,iAwake).GroupSpeed(iFOVneuro,:,iFun)=repmat(nanmean(AlignedSpeedMeta(:,I1),2)',length(iFOVneuro), 1);
% 
%                preSLMdata=squeeze(AlignedtempNData(:,(PSTHparam.PreSLMCal-PSTHparam.PreTestFrame+1):PSTHparam.PreSLMCal,I1));
%                postSLMdata=squeeze(AlignedtempNData(:,PSTHparam.PreSLMCal+[1:TestStepFrame],I1));
%                for jCell=1:sum(InputFOVData.NeuroPos3DMeta(:,4)==iFOV)
% 
%                     temp1=squeeze(nanmean(preSLMdata(jCell,:,:),2));
%                     temp2=squeeze(nanmean(postSLMdata(jCell,:,:),2));
%                     [~,p(jCell,1),~,t(jCell)]=ttest(temp2(:),temp1(:));
%                     temp3(jCell)=mean(temp2(:))-mean(temp1(:));
% 
%                end
%                statGroupRes(iVol,iAwake).p(iFOVneuro,iFun)=p;
%                % statGroupRes(iFun).t=t;
%                statGroupRes(iVol,iAwake).delta(iFOVneuro,iFun)=temp3;
%                statGroupRes(iVol,iAwake).GroupSampleN(iFOVneuro,iFun)=length(I1);
% 
%                clear p t temp3;
% 
% 
%             else
%                % statGroupRes(iFun).p=zeros(size(iscell))+1;
%                % statGroupRes(iFun).delta=zeros(size(iscell));
% 
%             end
% 
%         end
%  %%Add Zero power as addtional group
% 
%         statGroupRes(iVol,iAwake).deltaNorm=statGroupRes(iVol,iAwake).delta(:,length(GroupList))-repmat(statGroupRes(iVol,iAwake).delta(:,length(GroupList)+1),1,length(GroupList));
% 
%         statGroupRes(iVol,iAwake).GroupSpeedTrial=SpeedGroup;
%         clear SpeedGroup
% 
%     end
% end



