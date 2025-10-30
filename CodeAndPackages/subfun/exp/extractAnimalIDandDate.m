function [AnimalID, DateInfo] = extractAnimalIDandDate(ProcessFolder)
    % Function to extract AnimalID (SLxxxx or Lxxxxx) and DateInfo (MMDDYYYY) 
    % from ProcessFolder path

    % Regular expression to match either 'SL' followed by 4 digits or 'L' followed by 5 digits
    animalIDPattern = '(SL\d{4}|L\d{5})';

    % Regular expression to match the 8-digit date (MMDDYYYY)
    datePattern = '\d{8}';

    % Extract AnimalID using regular expression
    animalIDMatch = regexp(ProcessFolder, animalIDPattern, 'match');

    % Extract date using regular expression
    dateMatch = regexp(ProcessFolder, datePattern, 'match');

    % Check if both patterns were found
    if ~isempty(animalIDMatch) && ~isempty(dateMatch)
        AnimalID = animalIDMatch{1};
        DateInfo = dateMatch{1};
    else
        AnimalID = '';
        DateInfo = '';
        warning('Could not extract AnimalID or DateInfo from the provided ProcessFolder string.');
    end
end
