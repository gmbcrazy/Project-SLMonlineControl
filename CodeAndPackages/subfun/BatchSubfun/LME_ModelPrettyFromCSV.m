function modelStr = LME_ModelPrettyFromCSV(csvPath, spec, responseName, randEffectName)
% LME_ModelPrettyFromCSV  Build a compact LME formula like
%
%   Response ~ (A + B + C) * (D + E) + H + (1 | Session)
%
% from rows in an LME summary CSV, WITHOUT explicitly listing interaction
% terms (they are implied by the "*" operator).
%
% INPUTS
%   csvPath        : path to LME_FixedEffectsSummary_*.csv
%
%   spec           : struct with identifying fields; ONLY these are used:
%                      .DataType
%                      .Window
%                      .Condition
%                      .Subset
%                      .Model
%
%   responseName   : (optional) name of response; default 'Response'
%   randEffectName : (optional) grouping factor for random intercept;
%                    default 'Session'. If '' or [], no random part added.
%
% OUTPUT
%   modelStr       : char, compact formula string.
%

    if nargin < 3 || isempty(responseName)
        responseName = 'Response';
    end
    if nargin < 4
        randEffectName = 'Session';
    end

    T = readtable(csvPath);

    if ~ismember('Term', T.Properties.VariableNames)
        error('CSV must contain a column named "Term".');
    end

    % ------------------ filter rows by spec (model identity) ------------------ %
    mask = true(height(T),1);
    filterFields = {'DataType','Window','Condition','Subset','Model'};

    for i = 1:numel(filterFields)
        fn = filterFields{i};
        if isfield(spec, fn) && ismember(fn, T.Properties.VariableNames)
            val = spec.(fn);
            if ~isempty(val)
                if iscell(val), val = val{1}; end
                mask = mask & strcmp(T.(fn), val);
            end
        end
    end

    subT = T(mask,:);
    if isempty(subT)
        error('No rows in CSV match the given spec (DataType/Window/Condition/Subset/Model).');
    end

    % ------------------ collect all terms ------------------ %
    allTerms = unique(subT.Term, 'stable');

    % separate main vs interaction
    isInter = contains(allTerms, ':');
    mains   = allTerms(~isInter);
    inter   = allTerms(isInter);

    % remove intercept from main effects, if present
    mains = mains(~strcmp(mains, '(Intercept)'));

    % if no interactions, just return additive model
    if isempty(inter)
        rhs = strjoin(mains, ' + ');
        if ~isempty(randEffectName)
            rhs = sprintf('%s + (1 | %s)', rhs, randEffectName);
        end
        modelStr = sprintf('%s ~ %s', responseName, rhs);
        return;
    end

    % ------------------ decompose interactions into pairs ------------------ %
    pairs = cellfun(@(s) strsplit(s, ':'), inter, 'UniformOutput', false);
    pairs = cellfun(@(p) strtrim(p), pairs, 'UniformOutput', false);

    leftTokens  = cellfun(@(p) p{1}, pairs, 'UniformOutput', false);
    rightTokens = cellfun(@(p) p{2}, pairs, 'UniformOutput', false);

    % initial group guesses from positions
    G1 = unique(leftTokens,  'stable');
    G2 = unique(rightTokens, 'stable');

    % ------------------ classify interactions into main block vs extras ------------------ %
    inMainBlock = false(numel(inter),1); %#ok<NASGU>

    for i = 1:numel(inter)
        A = pairs{i}{1};
        B = pairs{i}{2};
        % interaction is considered part of the main block if
        % (A in G1 & B in G2) OR (A in G2 & B in G1)
        if (ismember(A,G1) && ismember(B,G2)) || ...
           (ismember(A,G2) && ismember(B,G1))
            inMainBlock(i) = true;
        else
            % previously we collected "extraInterStrings" here, but per request
            % we now IGNORE explicit interaction terms in the printed model.
        end
    end

    % ------------------ main effects covered by the (G1)*(G2) block ------------------ %
    mainsCovered   = unique([G1; G2],'stable');
    extraMainTerms = mains(~ismember(mains, mainsCovered));

    % ------------------ build compact RHS ------------------ %
    rhsParts = {};

    % main (G1)*(G2) block:
    group1Str = strjoin(G1, ' + ');
    group2Str = strjoin(G2, ' + ');
    rhsParts{end+1} = sprintf('(%s) * (%s)', group1Str, group2Str);

    % extra main effects:
    if ~isempty(extraMainTerms)
        rhsParts{end+1} = strjoin(extraMainTerms, ' + ');
    end

    % NO explicit "+ F * G" interaction terms are added here.

    % random part:
    if ~isempty(randEffectName)
        rhsParts{end+1} = sprintf('(1 | %s)', randEffectName);
    end

    rhs = strjoin(rhsParts, ' + ');

    % final formula
    modelStr = sprintf('%s ~ %s', responseName, rhs);
end
