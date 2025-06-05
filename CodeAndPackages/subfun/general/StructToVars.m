function StructToVars(s)
% Given a struct s, this creates a variable for each field in the caller workspace
fields = fieldnames(s);
for i = 1:length(fields)
    assignin('caller', fields{i}, s.(fields{i}));
end
end

% Usage:
%   StructToVars(ProcessPar);
% Now all fields of ProcessPar (e.g., GroupLabel, GroupColor, etc) are variables in your workspace.
