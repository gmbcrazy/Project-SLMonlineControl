function B = convertTableEntries(A)
    % Initialize an empty table B with the same size and variable names as A
    B = array2table(cell(size(A)), 'VariableNames', A.Properties.VariableNames);
    C=table2cell(A);
    A=cell2table(C,'VariableNames',A.Properties.VariableNames);
    % Loop through each column in A
    for col = 1:width(A)
        % Loop through each row in A
        for row = 1:height(A)
            % Extract the original string value from A
            originalValue = A.(col){row};
            if isempty(originalValue)
               B.(col){row}=[];
               continue
            end
            % Check if the original value is a string of numbers
            if isNumericString(originalValue)
                % Convert the string of numbers to a number
                B.(col){row} = str2double(originalValue);
            else
                % Keep the non-number string as is
                B.(col){row} = originalValue;
            end
        end
    end
    B=cell2table(table2array(B),'VariableNames',B.Properties.VariableNames);

end

function isNumeric = isNumericString(inputString)
    % Remove leading and trailing whitespaces
    inputString = strtrim(inputString);

    % Check if the string is not empty
    if isempty(inputString)
        isNumeric = false;
        return;
    end

    % Check if the string is a valid numeric string
    isNumeric = ~isnan(str2double(inputString));
end
