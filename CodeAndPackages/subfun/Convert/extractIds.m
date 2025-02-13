function ids = extractIds(folderNames)
    % Initialize output array
    ids = zeros(size(folderNames));
    
    for i = 1:numel(folderNames)
        % Extract the last three digits from the folder name
        tokens = regexp(folderNames{i}, '(\d{3})(?=\\?$)', 'tokens');
        
        if ~isempty(tokens)
            ids(i) = str2double(tokens{1}{1});
        else
            ids(i) = NaN; % Assign NaN if no match is found
        end
    end
end