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
