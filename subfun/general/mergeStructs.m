function C = mergeStructs(A, B)
    % Initialize structure C with the fields from A
    C = A;

    % Loop through the fields of B
    fieldsB = fieldnames(B);
    for i = 1:numel(fieldsB)
        % Check if the field already exists in C
        if isfield(C, fieldsB{i})
            % Field exists, concatenate the values
            C.(fieldsB{i}) = [C.(fieldsB{i}); B.(fieldsB{i})];
        else
            % Field doesn't exist, add the field and its value
            C.(fieldsB{i}) = B.(fieldsB{i});
        end
    end
end
