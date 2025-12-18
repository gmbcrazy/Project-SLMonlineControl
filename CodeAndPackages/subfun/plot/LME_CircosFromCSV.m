function [GgraphOut, theta] = LME_CircosFromCSV(csvPath, spec, plotCol, pThresh)
% LME_CircosFromCSV  Make a circos-style graph from LME summary CSV.
%
%   [GgraphOut, theta] = LME_CircosFromCSV(csvPath, spec, plotCol, pThresh)
%
%   Inputs
%   ------
%   csvPath : char, path to LME_FixedEffectsSummary_*.csv
%
%   spec    : struct describing which subset/model to plot. Typical fields:
%             - DataType    (e.g. 'deltaF')
%             - Window      (e.g. 'Win1')
%             - CellSubset  (e.g. 'All' or 'NonTarget')
%             - Condition   (e.g. 'NonSensory', 'Sensory')
%             - Subset      (e.g. 'NonSensory_G1')
%             - Model       (e.g. 'Base', 'Extended')
%             - ResponseVar (e.g. 'SpeedScore', 'SensoryScore')
%             - Term        (optional) cellstr/string array of term names
%                            to use as nodes (and their pairwise
%                            interactions)
%             - Clim        (optional) [min max] for colormap limits
%             - Colormap    (required) N×3 colormap used for edges/nodes
%             - NodeColor   (optional) default outline color for nodes
%             - TargetColor (optional) outline color for nodes whose
%                            name contains 'Target'
%
%   plotCol : char, numeric column to use as weight/color
%             (e.g. 'Estimate', 'tValue', etc.)
%
%   pThresh : scalar, p-value threshold. Only interactions with
%             pValue <= pThresh are plotted as edges. Main-effect
%             nodes with pValue <= pThresh are drawn bigger.
%
%   Outputs
%   -------
%   GgraphOut : struct with fields .G and .p (graph and plot handle)
%   theta     : angular position of nodes (from ClockWiseGraph)
%
%   Requires your functions:
%   - MarkovState_HeatStrPlot
%   - ClockWiseGraph
%

    % --------- read table ---------
    T = readtable(csvPath);

    if ~ismember('Term', T.Properties.VariableNames)
        error('CSV file must contain a column named ''Term''.');
    end
    if ~ismember('pValue', T.Properties.VariableNames)
        error('CSV file must contain a column named ''pValue''.');
    end
    if ~ismember(plotCol, T.Properties.VariableNames)
        error('Column "%s" not found in CSV.', plotCol);
    end
    if ~isnumeric(T.(plotCol))
        error('Column "%s" must be numeric.', plotCol);
    end

    % --------- filter rows for this model spec (except style fields) ---------
    mask = true(height(T),1);
    specFields = fieldnames(spec);
    
    for i = 1:numel(specFields)
        fn = specFields{i};
    
        % fields NOT used for filtering rows
        if ismember(fn, {'Term','Clim','Colormap','NodeColor','TargetColor'})
            continue;
        end
    
        % Only filter if BOTH:
        %   1) the CSV contains this column
        %   2) spec provides a non-empty value
        if ismember(fn, T.Properties.VariableNames) && ~isempty(spec.(fn))
            val = spec.(fn);
            if iscell(val)
                val = val{1};
            end
            mask = mask & strcmp(T.(fn), val);
        end
    end
    
    subT = T(mask,:);
    if isempty(subT)
        error(['No rows match the given spec (e.g. DataType/Window/CellSubset/' ...
               'Condition/Subset/Model/ResponseVar) in the CSV.']);
    end


    % --------- handle optional Term selection ---------
    if isfield(spec,'Term') && ~isempty(spec.Term)
        termList = spec.Term;
        if isstring(termList) || ischar(termList)
            termList = cellstr(termList);
        end
        termList = termList(:);
    else
        termList = [];  % means "use all terms found"
    end

    % --------- separate main effects vs interactions ---------
    isInteraction = contains(subT.Term, ':');

    % --- filter main effects by termList (if provided) ---
    if isempty(termList)
        mainT = subT(~isInteraction,:);
    else
        mainMask = ~isInteraction & ismember(subT.Term, termList);
        mainT    = subT(mainMask,:);
    end

    % --- filter interactions by termList (if provided) ---
    interT = subT(isInteraction,:);
    if ~isempty(termList)
        keepInt = false(height(interT),1);
        for i = 1:height(interT)
            parts = strsplit(interT.Term{i}, ':');
            parts = strtrim(parts);
            % keep interaction if *both* factors are in termList
            if all(ismember(parts, termList))
                keepInt(i) = true;
            end
        end
        interT = interT(keepInt,:);
    end

    % --------- collect node names ---------
    if isempty(termList)
        % derive from table
        nodeNames = mainT.Term;
        for i = 1:height(interT)
            parts = strsplit(interT.Term{i}, ':');
            nodeNames = [nodeNames; parts(:)];
        end
        nodeNames = unique(nodeNames, 'stable');
    else
        % use user-provided list (even if some are missing in mainT)
        nodeNames = termList;
    end
    nNodes = numel(nodeNames);

    % map name -> index
    name2idx = containers.Map(nodeNames, 1:nNodes);

    % --------- node values & p-values (from main effects if available) ---------
    nodeVal = zeros(nNodes,1);
    nodeP   = ones(nNodes,1)*NaN;

    for i = 1:nNodes
        thisName = nodeNames{i};
        idx = strcmp(mainT.Term, thisName);
        if any(idx)
            rows = mainT(idx,:);
            nodeVal(i) = mean(rows.(plotCol), 'omitnan');
            nodeP(i)   = mean(rows.pValue,    'omitnan');
        else
            nodeVal(i) = 0;
            nodeP(i)   = NaN;
        end
    end

    % --------- build adjacency matrix from significant interactions ---------
    Adj = zeros(nNodes);  % signed weight (for coloring)
    for i = 1:height(interT)
        termStr = interT.Term{i};
        parts = strsplit(termStr, ':');
        if numel(parts) ~= 2
            continue; % skip higher-order interactions if any
        end
        a = strtrim(parts{1});
        b = strtrim(parts{2});
        if ~isKey(name2idx,a) || ~isKey(name2idx,b)
            continue;
        end
        ia = name2idx(a);
        ib = name2idx(b);

        w  = interT.(plotCol)(i);
        pv = interT.pValue(i);

        if pv <= pThresh
            Adj(ia,ib) = w;
            % Adj(ib,ia) = w;  % treat as undirected if you want symmetry
        end
    end

    % --------- colormap limits: red = pos, blue = neg ---------
    if isfield(spec,'Clim') && ~isempty(spec.Clim)
        Clim = spec.Clim;
    else
        vals = [nodeVal(:); Adj(:)];
        vals = vals(~isnan(vals));
        if isempty(vals)
            vals = 0;
        end
        maxAbs = max(abs(vals));
        if maxAbs == 0
            maxAbs = 1;
        end
        Clim = [-maxAbs maxAbs];
    end
    ColorMapC = spec.Colormap;

    % --------- define structure & heat matrices for your circos ---------
    StrTransM  = double(Adj ~= 0);   % structure: edge exists or not
    HeatTransM = Adj;                % signed weights
    HeatTransM(StrTransM == 0) = NaN;  % ignore non-edges in color mapping

    HeatCounts = nodeVal;                 % node "heats"
    StrCounts  = zeros(nNodes,1) + 0.1;   % dummy positive
    ComColor   = ones(nNodes,3)*0.7;      % base node color

    % --------- call your circos plotter ---------
    [G,p] = MarkovState_HeatStrPlot(StrTransM,StrCounts, ...
                                    HeatTransM,HeatCounts, ...
                                    ComColor,ColorMapC,Clim);

    % ---- label nodes by term name instead of 1,2,3,... ----
    try
        p.NodeLabel = nodeNames;
    catch
        warning('Unable to set node labels.');
    end

    radius = 10;

    theta = ClockWiseGraph(G,p,radius);

    % ---- place labels *outside* the circle border ----
    % remove built-in labels so they don't sit on the nodes
    p.NodeLabel = repmat({''}, size(p.XData));

    hold on;
    labelRadius = radius * 1.50;   % outside the node circle
    set(gca,'xlim',[-1 1]*radius*2,'ylim',[-1 1]*radius*2);

    for i = 1:nNodes
        lx = labelRadius * cos(theta(i));
        ly = labelRadius * sin(theta(i));
        text(lx, ly, nodeNames{i}, ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','middle', ...
            'FontSize', 6);
    end

    % --------- adjust node sizes by significance ---------
    bigSize   = 8;    % significant node size
    smallSize = 0.1;  % non-significant node size

    mkSize = p.MarkerSize;  % current sizes (scalar or vector)
    if numel(mkSize) == 1
        mkSize = repmat(mkSize, nNodes,1);
    end
    for i = 1:nNodes
        if ~isnan(nodeP(i)) && nodeP(i) <= pThresh
            mkSize(i) = bigSize;
        else
            mkSize(i) = smallSize;
        end
    end
    p.MarkerSize = mkSize;

    % --------- make non-significant nodes white ---------
    nodeColor = p.NodeColor;
    for i = 1:nNodes
        if isnan(nodeP(i)) || nodeP(i) > pThresh
            nodeColor(i,:) = [1 1 1];
        end
    end
    p.NodeColor = nodeColor;

    % --------- output struct ---------
    GgraphOut.G = G;
    GgraphOut.p = p;
    GgraphOut.p.ArrowSize = 0;
    GgraphOut.p.LineWidth = ones(size(GgraphOut.p.LineWidth))*4;
end
