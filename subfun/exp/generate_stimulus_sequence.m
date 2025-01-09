function [sequence,CountStim] = generate_stimulus_sequence(StimNum, StimID, RepStimNumTh,MaxDiffNum)
    % Pre-allocate the sequence array
    sequence = zeros(1, StimNum);

    % Generate the sequence within constraints
    MakeSeq = 1;  % Control flag to continue generating until valid sequence is found

    while MakeSeq
        % Generate random sequence from StimID pool
        SequenceTemp = StimID(randi(length(StimID), StimNum, 1));
        CountStim=zeros(1,length(StimID));

        % Calculate differences between consecutive elements
        DiffSeq = diff(SequenceTemp);

        % Mark where the same stimulus occurs consecutively
        MarkSameSeq = (DiffSeq == 0);

        % Find periods of consecutive same stimuli
        SameSeqPeriod = MarkToPeriod(MarkSameSeq);
        SameSeqL = SameSeqPeriod(2, :) - SameSeqPeriod(1, :) + 2;

        % Check if the length of any same-stimulus period exceeds the threshold
        if max(SameSeqL) < RepStimNumTh
            if ~isempty(MaxDiffNum)
                for j=1:length(StimID)
                    CountStim(j)=sum(SequenceTemp==StimID(j));
                end
                if max(CountStim)-min(CountStim)<= MaxDiffNum
                   MakeSeq = 0;  % Valid sequence found
                   sequence = SequenceTemp';  % Save the generated sequence
                   break
                end

            else
                 MakeSeq = 0;  % Valid sequence found
                 sequence = SequenceTemp';  % Save the generated sequence
                 for j=1:length(StimID)
                     CountStim(j)=sum(SequenceTemp==StimID(j));
                 end
                break

            end

        end

    end


