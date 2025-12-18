function modelStr = LME_ModelFromCSV(csvPath, spec, responseName, randEffectName)
% LME_ModelFromCSV  Build a simple LME formula from the summary CSV:
%
%   Response ~ Term1 + Term2 + ... + TermK + (1 | Session)
%
% INPUTS
%   csvPath        : path to LME_FixedEffectsSummary_*.csv
%
%   spec           : struct used to filter rows; only these fields matter:
%                      .DataType
%                      .Window
%                      .Condition
%                      .Subset
%                      .Model
%                    Any of those present & non-empty will be used.
%                    All other fields (Term, Colormap, NodeColor, etc.)
%                    are ignored.
%
%   responseName   : (optional) name of response; default 'Response'
%   randEffectName : (optional) grouping factor for random intercept;
%                    default 'Session'.  If '' or [], no random term is added.
%
% OUTPUT
%   modelStr       : char, e.g.
%       'Response ~ SpeedR_z + SensoryR_z + TargetCellN_z + ... + (1 | Session)'
%

    if nargin < 3 || isempty(responseName)
        responseName = 'Response';
    end
    if nargin < 4
        randEffectName = 'Session';
    end

    % ------------------ read table ------------------ %
    T = readtable(csvPath);

    if ~ismember('Term', T.Properties.VariableNames)
        error('CSV must contain a column named "Term".');
    end

    % ------------------ filter by spec ------------------ %
    mask = true(height(T),1);

    % Only these spec fields are used to filter
    filterFields = {'DataType','Window','Condition','Subset','Model'};

    for i = 1:numel(filterFields)
        fn = filterFields{i};
        if isfield(spec, fn) && ismember(fn, T.Properties.VariableNames)
            val = spec.(fn);
            if ~isempty(val)
                if iscell(val), val = val{1}; end   % scalar in cell
                mask = mask & strcmp(T.(fn), val);
            end
        end
    end

    subT = T(mask,:);
    if isempty(subT)
        error('No rows in CSV match the given spec (DataType/Window/Condition/Subset/Model).');
    end

    % ------------------ collect ALL terms ------------------ %
    termsAll = unique(subT.Term, 'stable');  % keep original ordering

    % build RHS
    rhs = strjoin(termsAll, ' + ');

    % add random term if requested
    if ~isempty(randEffectName)
        rhs = sprintf('%s + (1 | %s)', rhs, randEffectName);
    end

    % final formula
    modelStr = sprintf('%s ~ %s', responseName, rhs);
end
