function visualizeLMESummary(csvPath, modelSpecs, colors, pThresh)
% visualizeLMESummary  Visualize and compare LME summary results from CSV.
%
%   visualizeLMESummary(csvPath, modelSpecs, colors, pThresh)
%
%   Inputs
%   ------
%   csvPath    : char/string, path to the summary CSV file
%
%   modelSpecs : struct array, one element per model to plot.
%                Allowed fields (all optional except PlotCol):
%                  DataType  - string, value in CSV column 'DataType'
%                  Window    - string, value in CSV column 'Window'
%                  Condition - string, value in CSV column 'Condition'
%                  Subset    - string, value in CSV column 'Subset'
%                  Term      - string or cellstr; if empty/missing,
%                              all terms from the CSV subset are used.
%                  PlotCol   - string, numeric column name to plot
%                              (e.g. 'Estimate', 'StdError',
%                              'tValue', 'pValue')
%
%   colors     : n x 3 RGB matrix; n must equal numel(modelSpecs).
%
%   pThresh    : scalar p-value threshold. Points with pValue > pThresh
%                are plotted as '.' (size 5); others as '*' (size 8).
%
%   The function creates a single figure with all models overlaid.
%   X-axis: term index (with labels). Y-axis: chosen PlotCol values.

    % ---------- read CSV ----------
    T = readtable(csvPath);

    nModels = numel(modelSpecs);
    if size(colors,1) ~= nModels
        error('colors must have one row per modelSpec (n x 3 RGB).');
    end

    % sanity check for pValue column
    if ~ismember('pValue', T.Properties.VariableNames)
        error('CSV file must contain a column named ''pValue''.');
    end

    % ---------- figure setup ----------
    hold on;
    legStr = cell(nModels,1);
    allTermsForAxis = {};

    for k = 1:nModels
        spec = modelSpecs(k);

        % ----- build logical mask for this model -----
        rows = true(height(T),1);

        fieldNames = {'DataType','Window','Condition','Subset'};
        for f = 1:numel(fieldNames)
            fn = fieldNames{f};
            if isfield(spec, fn) && ~isempty(spec.(fn))
                if ~ismember(fn, T.Properties.VariableNames)
                    error('Column ''%s'' not found in CSV.', fn);
                end
                rows = rows & strcmp(T.(fn), spec.(fn));
            end
        end

        subT = T(rows,:);
        if isempty(subT)
            warning('No rows matched for model %d; skipping.', k);
            continue;
        end

        % ----- determine term list for this model -----
        if isfield(spec,'Term') && ~isempty(spec.Term)
            if ischar(spec.Term) || isstring(spec.Term)
                termList = cellstr(spec.Term);
            else
                termList = spec.Term;     % assume cellstr
            end
        else
            termList = unique(subT.Term, 'stable');
        end

        % Use terms from first non-empty model to define x-axis
        if isempty(allTermsForAxis)
            allTermsForAxis = termList;
        end

        % restrict to terms present in this model AND in axis terms
        [commonTerms, iaAxis, ~] = intersect(allTermsForAxis, termList, 'stable');
        if isempty(commonTerms)
            warning('No overlapping Terms for model %d; skipping.', k);
            continue;
        end

        m = numel(commonTerms);

        % ----- validate PlotCol -----
        if ~isfield(spec,'PlotCol') || isempty(spec.PlotCol)
            error('modelSpecs(%d).PlotCol must be specified.', k);
        end
        plotCol = spec.PlotCol;
        if ~ismember(plotCol, subT.Properties.VariableNames)
            error('Column ''%s'' not found in CSV.', plotCol);
        end
        if ~isnumeric(subT.(plotCol))
            error('Column ''%s'' must be numeric.', plotCol);
        end

        yvals = nan(1,m);
        pvals = nan(1,m);
        xvals = iaAxis;  % indices into global term axis

        for i = 1:m
            thisTerm = commonTerms{i};
            idx = strcmp(subT.Term, thisTerm);
            if ~any(idx), continue; end

            rowsTerm = subT(idx,:);

            % if multiple rows (e.g. multiple windows) remain,
            % take the mean as a simple summary.
            vals = rowsTerm.(plotCol);
            yvals(i) = mean(vals, 'omitnan');

            pv = rowsTerm.pValue;
            pvals(i) = mean(pv, 'omitnan');
        end

        % ----- plot each point with marker based on p-value -----
        for i = 1:m
            if isnan(yvals(i)), continue; end
            if pvals(i) <= pThresh
                marker = '*';
                ms = 8;
            else
                marker = '.';
                ms = 5;
            end
            plot(xvals(i), yvals(i), marker, ...
                 'Color', colors(k,:), 'MarkerSize', ms);
        end

        % legend label: use Subset if available, otherwise DataType/Window/Condition
        if isfield(spec,'Subset') && ~isempty(spec.Subset)
            legStr{k} = spec.Subset;
        else
            parts = {};
            for f = 1:numel(fieldNames)
                fn = fieldNames{f};
                if isfield(spec, fn) && ~isempty(spec.(fn))
                    parts{end+1} = spec.(fn); %#ok<AGROW>
                end
            end
            legStr{k} = strjoin(parts,'_');
        end
    end

    % ----- finalize axes -----
    if isempty(allTermsForAxis)
        warning('No data plotted.');
        return;
    end

    set(gca, 'XTick', 1:numel(allTermsForAxis), ...
             'XTickLabel', allTermsForAxis, ...
             'XTickLabelRotation', 45);
    xlabel('Term');
    % Assume all models use the same PlotCol; just grab the first non-empty
    firstPlotCol = '';
    for k = 1:nModels
        if isfield(modelSpecs(k),'PlotCol') && ~isempty(modelSpecs(k).PlotCol)
            firstPlotCol = modelSpecs(k).PlotCol;
            break;
        end
    end
    if ~isempty(firstPlotCol)
        ylabel(firstPlotCol);
    else
        ylabel('Value');
    end
    % legend(legStr, 'Interpreter','none', 'Location','best');
    % box on; grid on;
    title('LME Summary Comparison');

    hold off;
end
