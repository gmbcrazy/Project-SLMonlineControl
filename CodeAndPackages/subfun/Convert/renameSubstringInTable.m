function T = renameSubstringInTable(T, oldStr, newStr)
% renameSubstringInTable  Replace substring in table variable names.
%
%   T = renameSubstringInTable(T, oldStr, newStr)
%
%   INPUTS:
%       T      : MATLAB table
%       oldStr : substring to replace
%       newStr : replacement substring
%
%   OUTPUT:
%       T      : table with updated variable names
%
%   Example:
%       T = renameSubstringInTable(T, 'Whisker', 'Sensory');

    if ~ischar(oldStr) && ~isstring(oldStr)
        error('oldStr must be a char or string.');
    end
    if ~ischar(newStr) && ~isstring(newStr)
        error('newStr must be a char or string.');
    end

    % Current variable names
    oldNames = T.Properties.VariableNames;
    newNames = oldNames;

    % Replace substring
    for i = 1:numel(oldNames)
        newNames{i} = strrep(oldNames{i}, oldStr, newStr);
    end

    % Apply new names
    T.Properties.VariableNames = newNames;
end
