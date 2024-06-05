% Define ySliderCallback to adjust the y-axis view of the plot
function ySliderCallback(src, eventdata)
    val = round(get(src, 'Value')); % Get the slider value, ensure it's negative
    ylim([val-2 val+2]); % Adjust y-axis to show cells around the slider value
end
