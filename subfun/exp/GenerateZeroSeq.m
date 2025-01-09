function SeqZero = GenerateZeroSeq(TotalSequence, StimID, ZeroN)
    % Pre-allocate the sequence of zeros
    SeqZero = zeros(size(TotalSequence));

    % Iterate for the number of zeros to insert
    for iProcess = 1:ZeroN
        % Count occurrences of each StimID in the sequence
        CountStim = zeros(1, length(StimID));
        for j = 1:length(StimID)
            CountStim(j) = sum(TotalSequence == StimID(j)&SeqZero==0);
        end

        % Find the stimulus ID with the maximum count
        [~, DelInd] = max(CountStim);
        DelID = StimID(DelInd);

        % Find indices where the stimulus ID repeats consecutively
        ZeroInd = find((diff([0; TotalSequence]) == 0 & TotalSequence == DelID) == 1);

        % If there are no consecutive repeats, randomly pick an index with the same StimID
        if isempty(ZeroInd)
            ZeroInd = find(TotalSequence == DelID);
            ZeroInd = ZeroInd(randi(length(ZeroInd), 1));
        else
            % Randomly pick one of the consecutive repeat indices
            ZeroInd = ZeroInd(randi(length(ZeroInd), 1));
        end

        % Mark the selected index as zero in SeqZero
        SeqZero(ZeroInd) = 1;
        % Decrement the count of zeros to insert
        ZeroN = ZeroN - 1;
    end
end
