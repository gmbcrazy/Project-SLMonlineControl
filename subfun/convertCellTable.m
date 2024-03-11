function convertedTable = convertCellTable(inputTable)
    % Copy the input table to the output table
    convertedTable = inputTable;
    convertedTable = createEmptyTable(inputTable);

    C=table2cell(inputTable);
    inputTable=cell2table(C,'VariableNames',inputTable.Properties.VariableNames);


    % Convert each cell item to either a number, a string, or a logical variable
    for col = 1:width(inputTable)
        columnData = inputTable.(col);
        for i = 1:numel(columnData)
            cellValue = columnData{i};
            
            % Try to convert to a number
            numValue = str2double(cellValue);
            if ~isnan(numValue)
                convertedTable.(col)(i) = numValue;
            else
                % Try to convert to a logical variable
                if strcmpi(cellValue, 'True') || strcmpi(cellValue, 'False')
                    convertedTable.(col)(i) = strcmpi(cellValue, 'True');
                end
                % If neither a number nor 'True'/'False', keep it as a string
            end
        end
    end
end


function A = createEmptyTable(B)
    % Create a new table A with empty entries
    emptyData = zeros(size(B));
    A = array2table(emptyData, 'VariableNames', B.Properties.VariableNames);
end