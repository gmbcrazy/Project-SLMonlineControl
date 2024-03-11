function C = addRowsTables(A, B)
    % Identify shared and non-shared variables
    sharedVars = intersect(A.Properties.VariableNames, B.Properties.VariableNames);
    nonSharedVarsA = setdiff(A.Properties.VariableNames, sharedVars);
    nonSharedVarsB = setdiff(B.Properties.VariableNames, sharedVars);

    % Reorder variables in A and B to have the same order
    A = A(:, [sharedVars, nonSharedVarsB]);
    B = B(:, [sharedVars, nonSharedVarsA]);

    % Add missing non-shared variables to C with zeros as placeholders
    for var = nonSharedVarsB
        B.(var{1}) = zeros(height(B), 1);
    end

    for var = nonSharedVarsA
        A.(var{1}) = zeros(height(A), 1);
    end

    % Vertically concatenate tables A and B
    C = [A; B];
end
