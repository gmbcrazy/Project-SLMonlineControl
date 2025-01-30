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
        if iscell(columnData)
        for i = 1:numel(columnData)
            strValue = columnData{i};
            if isstrprop(strValue, 'digit')
                tableWithNumeric.(col)(i) = str2double(strValue);
            % else
            %     tableWithNumeric.(col)(i)=strValue;
            end
        end
        elseif ischar(columnData)
            strValue = columnData(i,:);
            if isstrprop(strValue, 'digit')
                tableWithNumeric.(col)(i) = str2double(strValue);
            % else
            %     tableWithNumeric.(col)(i)=strValue;
            end
        else

        end
    end

    % tableWithNumeric = inputTable;
end

function A = createEmptyTable(B)
    % Create a new table A with empty entries
    emptyData = zeros(size(B));
    A = array2table(emptyData, 'VariableNames', B.Properties.VariableNames);
end