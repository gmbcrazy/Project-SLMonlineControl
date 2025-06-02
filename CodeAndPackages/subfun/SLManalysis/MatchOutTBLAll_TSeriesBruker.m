function OutTBLtemp = MatchOutTBLAll_TSeriesBruker(OutTBLAll, TSeriesBrukerTBL)
    % Function to integrate TSeriesBrukerTBL data into OutTBLAll based on FileID matching
    % using 'Group' (SynMPFunGroup) and cumulative 'markCycle' (Reps)
    % 
    % Inputs:
    %   OutTBLAll - Table containing experiment information
    %   TSeriesBrukerTBL - Cell array of tables containing additional experimental data
    % 
    % Output:
    %   OutTBLAll - Updated table with matched data
    
    uniqueFileIDs = unique(OutTBLAll.FileID);
    OutTBLtemp=[];
    for i = 1:length(uniqueFileIDs)
        fileID = uniqueFileIDs(i);
        subTable = OutTBLAll(OutTBLAll.FileID == fileID, :);
        
        % Try to find matching TSeriesBrukerTBL entry
        for j = 1:length(TSeriesBrukerTBL)
            tsTable = TSeriesBrukerTBL{j};
            tsTable.TSeriesInd=repmat(j,size(tsTable,1),1);
            % Check if Group matches SynMPFunGroup
            % MarkI1=strfind(tsTable.SynMPFunGroup',subTable.Group(:)');
            if strfind(tsTable.SynMPFunGroup',subTable.Group(:)')
                
                % Compute cumulative Reps
                Cycles = cumsum(tsTable.Reps)+1;
                if length(Cycles) > 1
                    Cycles = Cycles(1:end-1); % Ignore last cumulative sum element
                end
                
                % Check if markCycle matches Cycles
                MarkI=strfind(Cycles(:)',subTable.markCycle(:)');
                if  ~isempty(MarkI)
                    
                    % Assign matching TSeriesBrukerTBL information
                    matchedIdxTemp = ismember(OutTBLAll.FileID, fileID);
                    MarkI=MarkI+1;
                    AddTable=tsTable(MarkI:MarkI+sum(matchedIdxTemp)-1,:);
                    if exist('matchedIdx')
                       matchedIdx=matchedIdx+matchedIdxTemp;
                    else
                       matchedIdx=matchedIdxTemp;
                    end

                    % ismember(subTable.markCycle, Cycles)

                    % OutTBLAll(matchedIdx, :) = [subTable, AddTable];
                    OutTBLtemp=[OutTBLtemp;subTable, AddTable];

% Now concatenate
                    % OutTBLtemp = [OutTBLtemp;subTable];
                    break;
                end
            end
        end

    end


end
