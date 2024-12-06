function numericTable = convertStringsToNumbers(inputTable)
    % Copy the input table to the output table

    % Convert strings like 'True' to 1 and 'False' to 0
    numericTable = convertLogicalStrings(inputTable);

    % Convert strings representing numeric values to numbers
    numericTable = convertNumericStrings(numericTable);
end


function A = createEmptyTable(B)
    % Create a new table A with empty entries
    emptyData = zeros(size(B));
    A = array2table(emptyData, 'VariableNames', B.Properties.VariableNames);
end

function tableWithLogical = convertLogicalStrings(inputTable)
    % Convert 'True' to 1 and 'False' to 0
    tableWithLogical=createEmptyTable(inputTable);
    logicalStrings = {'True', 'False'};
    for col = 1:width(inputTable)
        columnData = inputTable.(col);
        for i = 1:numel(columnData)
            if ismember(columnData{i}, logicalStrings)
                tableWithLogical.(col)(i) = find(strcmp(columnData{i}, logicalStrings)) - 1;
            % else
            %     tableWithLogical.(col)(i) = columnData{i};
                
            end
        end
    end

    % tableWithLogical = inputTable;
end

function tableWithNumeric = convertNumericStrings(inputTable)
    % Convert strings representing numeric values to numbers
        tableWithNumeric=createEmptyTable(inputTable);

    for col = 1:width(inputTable)
        columnData = inputTable.(col);
        if isnumeric(columnData)
           for i = 1:numel(columnData)
               tableWithNumeric.(col)(i)=columnData(i);
           end
           continue

        end
        for i = 1:numel(columnData)
            strValue = columnData{i};
            if isstrprop(strValue, 'digit') || (strValue(1) == '-' && isstrprop(strValue(2:end), 'digit'))
                tableWithNumeric.(col)(i) = str2double(strValue);
            % else
            %     tableWithNumeric.(col)(i)=strValue;
            end
        end
    end

    % tableWithNumeric = inputTable;
end
