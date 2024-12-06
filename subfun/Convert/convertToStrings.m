function B = convertToStrings(A)
    % Initialize the output cell variable
    B = cell(size(A));
    
    % Loop through each element of A
    for i = 1:numel(A)
        % Check if the element is already a string
        if ischar(A{i})
            % If it's a string, keep it as is
            B{i} = A{i};
        else
            % If it's not a string, convert it to a string
            B{i} = num2str(A{i});
        end
    end
end
