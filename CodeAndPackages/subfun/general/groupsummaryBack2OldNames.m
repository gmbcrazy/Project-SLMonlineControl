function newtbl2 = groupsummaryBack2OldNames(oldtbl, newtbl, GroupMethod)

    % Remove GroupCount column if present
    % if any(strcmp(newtbl.Properties.VariableNames, 'GroupCount'))
    %     newtbl(:, 'GroupCount') = [];   
    % end
    
    newtbl2 = newtbl;

    % Extract current variable names
    newNames = newtbl.Properties.VariableNames;

    % Find which names start with the GroupMethod prefix (e.g. 'mean_')
    processVars = startsWith(newNames, [GroupMethod '_']);

    % Remove the prefix only from those vars
    strippedNames = erase(newNames(processVars), [GroupMethod '_']);

    % Replace with stripped names
    newNames(processVars) = strippedNames;

    % Ensure uniqueness (avoid duplicate names like 'Session')
    newNames = matlab.lang.makeUniqueStrings(newNames);

    % Assign back
    newtbl2.Properties.VariableNames = newNames;
end
