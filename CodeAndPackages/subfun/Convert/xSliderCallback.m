% Define xSliderCallback to adjust the x-axis view of the plot
function xSliderCallback(src, eventdata)
    val = round(get(src, 'Value')); % Get the slider value
    xlim([max(1, val-500) val+500]); % Adjust x-axis to show around ±500 points from the slider value
end
