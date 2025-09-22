function [deltaFoF,deltaFoF2]=F2deltaFoF(F,Fneu,fR)

        F = F'; Fneu = Fneu';
        Fsubtracted = F - 0.7*Fneu;
        
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
        deltaFoF2   = bsxfun(@rdivide, deltaFoF , 1e-12 + sd); % THIS IS NOT deltaFoF; I JUST DID NOT WANT TO CHANGE ALL MY CODES
