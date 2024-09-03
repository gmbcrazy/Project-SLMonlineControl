function paddedStr = padNumber(num, totalLength)
    % padNumber - Converts a number to a string with leading zeros
    %
    % Syntax: paddedStr = padNumber(num, totalLength)
    %
    % Inputs:
    %    num - The number to be converted
    %    totalLength - The total length of the resulting string
    %
    % Outputs:
    %    paddedStr - The number converted to a string with leading zeros

    % Convert the number to a string with leading zeros
    paddedStr = sprintf(['%0' num2str(totalLength) 'd'], num);
end
