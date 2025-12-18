function stats=ErrorBarPlotPairTBL(xGroup, DataTBL, VariForCompare, GroupingVari, PlotColor, PlotMode, varargin)

%%%Input:
%%%%%%%xGroup: plotting order and x-axis coordinates, e.g., [2 1 3] means plot group 2, then 1, then 3
%%%%%%%DataTBL: table containing the data
%%%%%%%VariForCompare: column name in DataTBL for the variable to compare
%%%%%%%GroupingVari: column name in DataTBL for grouping (e.g., 'Group')
%%%%%%%PlotColor: color specification for groups
%%%%%%%PlotMode: 1=bar, 2=dot, 3=area, 4=dot+line, 5=line only
%%%%%%%varargin{1}: SessionVari - column name for session/subject ID (default: 'Session')
%%%%%%%varargin{2}: SavePath - path to save statistics (default: [])
%%%%%%%varargin{3}: GroupPair - structure with plotting parameters
%%%%%%%varargin{4}: AdjacentOnly - 1 to connect only adjacent pairs (default: 1)

% Parse optional inputs
if nargin >= 7
    SessionVari = varargin{1};
else
    SessionVari = 'Session';
end

if nargin >= 8
    SavePath = varargin{2};
else
    SavePath = [];
end

if nargin >= 9
    GroupPair = varargin{3};
else
    GroupPair = [];
end

if nargin >= 10
    AdjacentOnly = varargin{4};
else
    AdjacentOnly = 1; % Only connect adjacent pairs by default
end

% Initialize GroupPair with defaults
if isempty(GroupPair)
    GroupPair.CorrName = 'fdr';
    GroupPair.Q = 0.1;
    GroupPair.SignY = [];
    GroupPair.Plot = 1;
    GroupPair.Std = 1;
    GroupPair.SamplePlot = 1;
    GroupPair.SamplePairedPlot = 1;
    GroupPair.LimY = [];
    GroupPair.Marker = {'o'};
end

% Set up group names
uniqueGroups = unique(DataTBL.(GroupingVari));
nGroups = length(xGroup);

if ~isfield(GroupPair, 'GroupName')
    for i = 1:nGroups
        GroupPair.GroupName{i} = ['Group' num2str(xGroup(i))];
    end
end
GName = GroupPair.GroupName;

% Set up markers
if isfield(GroupPair, 'Marker')
    Marker = GroupPair.Marker;
else
    Marker = {'o'};
end
if length(Marker) == 1
    Marker = repmat(Marker, 1, nGroups);
end

% Handle color inputs
if isnumeric(PlotColor)
    barColor = PlotColor;
    barFaceColor = PlotColor;
elseif iscell(PlotColor)
    barColor = PlotColor{1};
    barFaceColor = PlotColor{2};
end

if size(barColor, 1) == 1
    barColor = repmat(barColor, nGroups, 1);
end
if size(barFaceColor, 1) == 1
    barFaceColor = repmat(barFaceColor, nGroups, 1);
end

% Extract data for each group
y_mean = cell(1, nGroups);
for i = 1:nGroups
    groupIdx = DataTBL.(GroupingVari) == xGroup(i);
    y_mean{i} = DataTBL.(VariForCompare)(groupIdx);
end

% Calculate mean and std/ste for each group
for ii = 1:nGroups
    tempMean(ii) = nanmean(y_mean{ii});
    if GroupPair.Std == 1
        tempStd(ii) = nanstd(y_mean{ii});
    else
        tempStd(ii) = nanstd(y_mean{ii}) / sqrt(sum(~isnan(y_mean{ii})));
    end
end

% Set up x-axis positions
x = 1:nGroups;
deltaX = 1;

% Organize data by session for paired plotting
sessionIDs = unique(DataTBL.(SessionVari));
nSessions = length(sessionIDs);

% Create matrices for paired data
xTemp = cell(1, nGroups);
yTemp = cell(1, nGroups);
sessionTemp = cell(1, nGroups);

for ii = 1:nGroups
    groupIdx = DataTBL.(GroupingVari) == xGroup(ii);
    Num(ii) = sum(groupIdx);
    yTemp{ii} = DataTBL.(VariForCompare)(groupIdx);
    sessionTemp{ii} = DataTBL.(SessionVari)(groupIdx);
    xTemp{ii} = randn(Num(ii), 1) * deltaX / 40 + x(ii);
end

% Plot individual samples
if GroupPair.SamplePlot == 1
    for ii = 1:nGroups
        plot(xTemp{ii}, yTemp{ii}, 'color', [0.6 0.6 0.6], 'linestyle', 'none', 'marker', '.', 'markersize', 5);
        hold on;
    end
end

% Plot connecting lines for paired samples
if GroupPair.SamplePlot == 1 && GroupPair.SamplePairedPlot ~= 0
    if AdjacentOnly == 1
        % Only connect adjacent groups
        pairsToConnect = [(1:nGroups-1)', (2:nGroups)'];
    else
        % Connect all pairs
        pairsToConnect = [];
        for i = 1:nGroups
            for j = (i+1):nGroups
                pairsToConnect = [pairsToConnect; i, j];
            end
        end
    end
    
    for iPair = 1:size(pairsToConnect, 1)
        ig1 = pairsToConnect(iPair, 1);
        ig2 = pairsToConnect(iPair, 2);
        
        % Find matching sessions
        [commonSessions, idx1, idx2] = intersect(sessionTemp{ig1}, sessionTemp{ig2});
        
        for iMatch = 1:length(commonSessions)
            x1 = xTemp{ig1}(idx1(iMatch));
            y1 = yTemp{ig1}(idx1(iMatch));
            x2 = xTemp{ig2}(idx2(iMatch));
            y2 = yTemp{ig2}(idx2(iMatch));
            
            if GroupPair.SamplePairedPlot == 1
                % Simple gray line
                plot([x1 x2], [y1 y2], 'color', [0.6 0.6 0.6], 'linestyle', '-', 'linewidth', 0.4);
            elseif GroupPair.SamplePairedPlot == 2
                % Red for increase, blue for decrease
                if y1 < y2
                    plot([x1 x2], [y1 y2], 'color', [1 0.5 0.5], 'linestyle', '-', 'linewidth', 0.4);
                elseif y1 > y2
                    plot([x1 x2], [y1 y2], 'color', [0.5 0.5 1], 'linestyle', '-', 'linewidth', 0.4);
                else
                    plot([x1 x2], [y1 y2], 'color', [0.6 0.6 0.6], 'linestyle', '-', 'linewidth', 0.4);
                end
            end
        end
    end
end

% Plot mean and error bars
if PlotMode == 1
    for i = 1:nGroups
        bar(x(i), tempMean(i), 0.8, 'facecolor', barFaceColor(i,:), 'Edgecolor', barColor(i,:)); hold on;
        plot([x(i) x(i)], tempMean(i) + [-tempStd(i) tempStd(i)], 'color', barColor(i,:));
    end
elseif PlotMode == 2
    for i = 1:nGroups
        plot(x(i), tempMean(i), 'color', barColor(i,:), 'linestyle', 'none', 'marker', Marker{i}, 'markersize', 5, 'markerfacecolor', barFaceColor(i,:)); hold on
        plot([x(i) x(i)], tempMean(i) + [-tempStd(i) tempStd(i)], 'color', barColor(i,:), 'linewidth', 1);
    end
elseif PlotMode == 4
    for i = 1:nGroups
        plot(x(i), tempMean(i), 'color', barColor(i,:), 'linestyle', 'none', 'marker', Marker{i}, 'markersize', 5, 'markerfacecolor', barFaceColor(i,:)); hold on
        plot([x(i) x(i)], tempMean(i) + [-tempStd(i) tempStd(i)], 'color', barColor(i,:), 'linewidth', 1);
    end
    plot(x, tempMean, 'color', barColor(end,:), 'linestyle', '-'); hold on
elseif PlotMode == 5
    plot(x, tempMean, 'color', barColor(end,:), 'linewidth', 1);
    plot(x, tempMean - tempStd, 'color', sqrt(barColor(end,:)), 'linewidth', 0.5);
    plot(x, tempMean + tempStd, 'color', sqrt(barColor(end,:)), 'linewidth', 0.5);
end

% Statistical testing - paired t-tests only
ShowTemp(1,:) = '                                                                        ';
ShowTemp(end+1, 1:length('Paired T-Tests')) = 'Paired T-Tests';

p = [];
Tstat = struct([]);
pairIndices = [];

% Determine which pairs to test
if AdjacentOnly == 1
    % Only test adjacent groups
    for i = 1:(nGroups-1)
        pairIndices = [pairIndices; i, i+1];
    end
else
    % Test all pairs
    for i = 1:nGroups
        for j = (i+1):nGroups
            pairIndices = [pairIndices; i, j];
        end
    end
end

% Perform paired t-tests
for iPair = 1:size(pairIndices, 1)
    ii = pairIndices(iPair, 1);
    iii = pairIndices(iPair, 2);
    
    % Find matching sessions
    [commonSessions, idx1, idx2] = intersect(sessionTemp{ii}, sessionTemp{iii});
    
    if length(commonSessions) > 1
        yTemp1 = yTemp{ii}(idx1);
        yTemp2 = yTemp{iii}(idx2);
        
        % Remove NaN pairs
        validIdx = ~isnan(yTemp1) & ~isnan(yTemp2);
        yTemp1 = yTemp1(validIdx);
        yTemp2 = yTemp2(validIdx);
        
        if length(yTemp1) > 1
            if iPair==1
            [~, p(iPair), ~, Tstat] = ttest(yTemp1, yTemp2);
            else
            [~, p(iPair), ~, Tstat(iPair)] = ttest(yTemp1, yTemp2);
            end
            TempText = [GName{ii} '-' GName{iii}, ' ,t' num2str(Tstat(iPair).df) '=' num2str(Tstat(iPair).tstat) ',p=' num2str(p(iPair)) ' (n=' num2str(length(yTemp1)) ' pairs)'];
        else
            p(iPair) = NaN;
            TempText = [GName{ii} '-' GName{iii}, ' - insufficient paired data'];
        end
    else
        p(iPair) = NaN;
        TempText = [GName{ii} '-' GName{iii}, ' - no paired data'];
    end
    ShowTemp(end+1, 1:length(TempText)) = TempText;
end

% FDR correction
validP = p(~isnan(p));
if ~isempty(validP)
    if strcmp(GroupPair.CorrName, 'fdr')
        if isfield(GroupPair, 'Crit_p')
            crit_p = GroupPair.Crit_p;
        else
            [h, crit_p, adj_p] = fdr_bh(validP, GroupPair.Q, 'pdep', 'yes');
        end
        crit_p(crit_p > 0.05) = 0.05;
        
        TempText = ['Multi-Comparison ' GroupPair.CorrName];
        ShowTemp(end+1, 1:length(TempText)) = TempText;
        TempText = ['Q=' num2str(GroupPair.Q) ' Pth=' num2str(crit_p)];
        ShowTemp(end+1, 1:length(TempText)) = TempText;
        
        % Plot significance markers
        if GroupPair.Plot == 1
            if isempty(GroupPair.SignY)
                GroupPair.SignY = max(tempMean) + max(tempStd) * 1.5;
            end
            if isempty(GroupPair.LimY)
                GroupPair.LimY = [min(tempMean) - max(tempStd) * 1.5, max(tempMean) + max(tempStd) * 1.5];
            end
            
            for iPair = 1:length(p)
                if ~isnan(p(iPair)) && p(iPair) <= crit_p
                    ii = pairIndices(iPair, 1);
                    iii = pairIndices(iPair, 2);
                    xP = [x(ii) + deltaX/20, x(iii) - deltaX/20];
                    yP = GroupPair.SignY + (iii - ii - 1) * abs(diff(GroupPair.LimY)) / 20;
                    
                    plot(xP, [yP yP], 'k-');
                    text(mean(xP), yP, '*', 'horizontalalignment', 'center', 'verticalalignment', 'baseline', 'fontsize', 10);
                end
            end
        end
    end
end

% Summary statistics
TempText = 'mean ± sd (n)';
ShowTemp(end+1, 1:length(TempText)) = TempText;
for ii = 1:nGroups
    TempText = [GName{ii} ' ' num2str(tempMean(ii)) ' ± ' num2str(nanstd(y_mean{ii})) ' (n=' num2str(sum(~isnan(y_mean{ii}))) ')'];
    ShowTemp(end+1, 1:length(TempText)) = TempText;
end

% Store statistics
stats.p_Pair = p;
stats.t_Pair = Tstat;
stats.pairIndices = pairIndices;
stats.GroupNames = GName;
stats.mean = tempMean;
stats.std = arrayfun(@(i) nanstd(y_mean{i}), 1:nGroups);
if ~isempty(validP)
    stats.crit_p = crit_p;
end

% Save to file
if ~isempty(SavePath)
    fileID = fopen(SavePath, 'w');
    for i = 1:size(ShowTemp, 1)
        fprintf(fileID, '%s\r\n', deblank(ShowTemp(i,:)));
    end
    fclose(fileID);
end

% Set axis properties
set(gca, 'XTick', x, 'XTickLabel', GName);
if ~isempty(GroupPair.LimY)
    ylim(GroupPair.LimY);
end