%%
% --- Ana revision 05/18/2022 --- removed excelSheetNum from input vars
% --- Ana revision 05/12/2022 --- xlsread replacement
% --- Ana revision 03/16/2022 --- Added spks norm schemes
% --- Ana revision 09/21/2021 ---

function[] = Merge(processedDir, processedDirSuite2P, excelfilename, rowsOfInterest)

% --- Set paths ---
expFiles = dir(processedDir);
expFiles = expFiles(~strncmpi('.', {expFiles.name}, 1));
expFilesSuite2P = dir(processedDirSuite2P);
expFilesSuite2P = expFilesSuite2P(~strncmpi('.', {expFilesSuite2P.name}, 1));

% [~,text,~] = xlsread(excelfilename,excelSheetNum);
opts = detectImportOptions(excelfilename);
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

clear i opts header text

for m = 1:numel(rowsOfInterest)
    
    k = rowsOfInterest(m)-shift;
    
    % --- Get Ana's data ---
    currExp = expFiles(contains({expFiles.name}, param.Animal{k}) == 1 & contains({expFiles.name}, param.Date{k}) == 1 & contains({expFiles.name}, param.ExpDivision{k}) == 1);
    numPlanes = size(currExp,1);
    
    for plane = 1:numPlanes
        
        % --- Get Ana's data ---
        currExpPlane = currExp(plane);
        currExpVars = dir([currExpPlane.folder '\' currExpPlane.name '\*vars.mat']);
        load([currExpVars.folder '\' currExpVars.name], 'caTrials', 'lastFrame', 'timeStampCa')
        
        fRTemp = [];
        for i = 1:length(caTrials)
            fRTemp(i) = mean(diff(timeStampCa(lastFrame(i)+1:lastFrame(i+1))));
        end
        fR = 1/mean(fRTemp); clear i fRTemp % fR in Hz
        
        % --- Get Suite2P data ---
        currExpPlaneSuite2P = expFilesSuite2P(contains({expFilesSuite2P.name}, param.Animal{k}) == 1 & contains({expFilesSuite2P.name}, param.Date{k}) == 1);
        currExpVarsSuite2P = dir([currExpPlaneSuite2P.folder '\' currExpPlaneSuite2P.name '\suite2p\plane' num2str(plane-1) '\Fall.mat']);
        iscell = [];
        load([currExpVarsSuite2P.folder '\' currExpVarsSuite2P.name])

        % --- Extract Suite2P data ---
        % On a second thought, it might be better just to keep everything
        % Just invert matrices because on my codes/Francisco's take m'
        % cellInd = find(iscell==1);
        % numberOfCells = length(cellInd);
        % redInd = redcell(cellInd,:);
        % Fcell = F(cellInd,:)';
        % Fneurop = Fneu(cellInd,:)';
        % Fdeconv = spks(cellInd,:)';
        F = F'; Fneu = Fneu'; spks = spks'; 
        Fsubtracted = F - 0.7*Fneu;
        stat=stat';
        ops=ops';
        % --- Calculate Baseline ---
        % Option 1 - fixed baseline (not ideal for complexe and long
        % experiments)
        % percentileCutOff = 10;
        % for c = 1:numberOfCells
        %    F0(:,c) = prctile(Fsubtracted(:,c),percentileCutOff); % I can use 3rd input as dimention to avoid cycle
        % end
        
        % Option 2 - running baseline
        % Hard to make any running baseline work well for high firing cells
        % which leads to option 3
        % percentileCutOff = 50;
        % win = round(fR*120); % 120s
        % F0 = zeros(size(Fsubtracted));
        % for c = 1:numberOfCells
        %    F0(:,c) = running_percentile(Fsubtracted(:,c), win, percentileCutOff);
        % end
        
        % lowPassFiltCutOff = 60; % in Hz
        % for c = 1:numberOfCells
        %    F0(:,c) = baselinePercentileFilter(Fsubtracted(:,c),fR,lowPassFiltCutOff,percentileCutOff);
        %end
        
        % Option 3 - remove peaks before estimating baseline, can do it
        % iteratively (see https://pubmed.ncbi.nlm.nih.gov/28132825/);
        % did not finish
        % comp = zeros(size(Fsubtracted));
        % for c = 1:numberOfCells
        %    d = [];
        %    d = [0; diff(Fsubtracted(:,c))];
        %    threshold = 3*minStd(d,round(fR*60));
        %    peaks1 = LocalMinima(d,1,-threshold);
        %    peaks2 = LocalMinima(-d,1,-threshold);
        %    comp(:,c) = 0;
        %    comp(peaks1,c) = 1;
        %    comp(peaks2,c) = 1;
        % end
        % conc = sort([peaks1; peaks2]);
        % dt = [0; diff(conc)];
        % shifted = dt <= 30;
        % shifted = [false;  shifted(2:end)];
        % tags = diff([shifted; false]);
        % starts = tags>0;
        % ends = tags<0;
        % onsets = conc(starts);
        % offsets = conc(ends);
        
        % Option 4 - from suite2p
        [NT , ~] = size(Fsubtracted);
        
        % determine and subtract the running baseline
        win         = 60; % in sec
        Fbase       = zeros(size(Fsubtracted), 'single');
        
        % if getOr(ops, 'runningBaseline', 0)
        Fbase   = Fsubtracted;
        %    if getOr(ops, 'zBaseline', 0)
        %        [~, ix] = sort(ops.zdrift);
        %        Fbase = Fbase(ix, :);
        %    end
        
        ntBase  = 2*ceil(win * fR/2)+1;
        Fbase   = cat(1, Fbase((ntBase-1)/2:-1:1, :), Fbase, Fbase(end:-1:end-(ntBase-1)/2, :));
        
        Fbase   = my_conv2(Fbase, 20, 1);
        Fbase   = movmin(Fbase, ntBase,1);
        Fbase   = movmax(Fbase, ntBase,1);
        
        Fbase   = Fbase((ntBase-1)/2 + [1:NT], :);
        
        %   if getOr(ops, 'zBaseline', 0)
        %       Fbase(ix,:) = Fbase;
        %   end
        % end
        
        deltaFoF    = Fsubtracted - Fbase; % THIS IS NOT deltaFoF; I JUST DID NOT WANT TO CHANGE ALL MY CODES
        
        % Fsort     = my_conv2(F1, ceil(ops.fs), 1);
        % Fsort     = sort(Fsort, 1, 'ascend');
        % baselines = Fsort(ceil(NT/20), :);
        % Fbase2    = bsxfun(@times, ones(NT,1), baselines);
        % F1        = F1 - Fbase2;
        % Fbase     = Fbase + Fbase2;
        
        % F1        = F1./max(std(F1,1,1), Fbase);
        
        % --- Normalization ---
        % Option 1 - DeltaF/F (can lead to a lot of distortions)
        % deltaFoF = zeros(size(Fsubtracted));
        % for c = 1:numberOfCells
        %   deltaFoF(:,c) = (Fsubtracted(:,c)-(F0(:,c)))./abs(F0(:,c));
        % end
        %
        % Norm - option 2 - From suite2p
        % normalize signal
        sd         = 1/sqrt(2) * std(deltaFoF(2:end, :) - deltaFoF(1:end-1, :), [], 1);
        deltaFoF   = bsxfun(@rdivide, deltaFoF , 1e-12 + sd); % THIS IS NOT deltaFoF; I JUST DID NOT WANT TO CHANGE ALL MY CODES

        % --- Append vars ---
        save([currExpVars.folder '\' currExpVars.name], 'iscell', 'redcell', 'F', 'Fneu', 'Fbase', 'deltaFoF', 'spks', 'fR','stat','ops', '-append');
    end
end