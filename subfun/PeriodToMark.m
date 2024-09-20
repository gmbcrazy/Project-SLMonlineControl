function Period = PeriodToMark(Period, MarkLength)
    % PerodToMark - Convert periods of occurrences to a binary marker vector
    %
    % This function takes in a matrix of periods and a length for the marker
    % vector, and converts the periods to a binary marker vector where 1 indicates
    % an occurrence and 0 indicates no occurrence.
    %
    % Inputs:
    %   Period - A 2-row matrix where Period(1, :) contains the start indices 
    %            of each period and Period(2, :) contains the end indices of each period.
    %   MarkLength - The length of the output marker vector.
    %
    % Output:
    %   Period - The input Period matrix is returned unchanged.
    
    % Initialize the marker vector with zeros
    Mark = zeros(1, MarkLength);
    
    % Loop through each period and set the corresponding indices in the marker vector to 1
    for iP = 1:size(Period, 2)
        Mark(Period(1, iP):Period(2, iP)) = 1;
    end
    
    %%%%%%% Mark is a vector with 0 (non-occurrence) and 1 (occurrence) values
    %%%%%%% Period gives the occurrence periods
    %%%%%%% Period(1, :) is the start index of each period
    %%%%%%% Period(2, :) is the end index of each period
    
    % Output the updated marker vector
    Period = Mark;
end
