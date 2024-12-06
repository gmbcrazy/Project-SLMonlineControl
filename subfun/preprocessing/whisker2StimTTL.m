function TTL = whisker2StimTTL(whiskerData, XTimesStdTh, MinInterVal)
    % whisker2StimTTL - Convert whisker data to TTL (transistor-transistor logic) signals for stimulation
    %
    % This function processes whisker data to generate TTL signals based on the data exceeding a specified threshold.
    %
    % Inputs:
    %   whiskerData - Matrix containing the whisker data.
    %   XTimesStdTh - Multiplier for the standard deviation to set the threshold.
    %   MinInterVal - Minimum interval between periods for merging.
    %
    % Output:
    %   TTL - TTL signal vector generated from the processed whisker data.

    % Process and normalize the whisker data
    temp = whiskerData(:, end); % Extract the last column of whiskerData
    temp = smooth2005((temp - min(temp)), 20, 'lowess'); % Smooth and normalize the data
    temp = AmpNormalize(temp); % Further normalize the amplitude
    
    % Plot threshold line
    % figure; % Create a new figure for plotting
    % plot([1 length(temp)], repmat(mean(temp) + XTimesStdTh * std(temp), 1, 2))
    
    % Calculate the threshold and create a marker vector based on the threshold
    threshold = mean(temp) + XTimesStdTh * std(temp);
    MarkP = temp > threshold; % Marker vector indicating where data exceeds threshold
    
    % Convert marker vector to periods and merge periods based on the minimum interval
    Period = MarkToPeriod(MarkP);
    Period2 = MergePeriod(Period, MinInterVal);
    Period3 = Period2;
    Period3(2,:)=Period3(1,:)+8;
    
    % Convert the merged periods back to a marker vector
    MarkP2 = PeriodToMark(Period3, length(MarkP));
    
    % Plot the processed data and marker vector
    % plot(temp); hold on;
    % plot(MarkP2);
    
    % Assign the final marker vector to the output TTL signal
    TTL = MarkP2;
end
