function AddCircosBand(CircosBand)

hold on;
for i=1:length(CircosBand)
    ColorOut = Value2Color(CircosBand(i).Values, CircosBand(i).Colormap, CircosBand(i).Clim);
    AddCircosBandSub(CircosBand(i).theta,CircosBand(i).rband,CircosBand(i).arc_gap,ColorOut);
    drawnow;
    set(gca, 'SortMethod', 'depth')

end

function AddCircosBandSub(theta,rband,arc_gap,rcolor)

r_outer = rband(1); % Outer radius (for heatmap band)
r_band_outer=rband(2);
N = length(theta);
% Define heatmap values (random example)
heatmap_values = rand(1, N); % Assign heat values (e.g., connectivity strength)

% Convert polar coordinates to Cartesian
x_outer = r_outer * cos(theta);
y_outer = r_outer * sin(theta);

hold on;
r_band_outer=r_outer+1;
% figure;
arc_width = (2*pi)/N;
for k = 1:N
    theta_arc = linspace(theta(k)-arc_width/2+arc_gap/2, theta(k)+arc_width/2-arc_gap/2, 30);
    x_arc_inner = r_outer * cos(theta_arc);
    y_arc_inner = r_outer * sin(theta_arc);
    x_arc_outer = r_band_outer * cos(fliplr(theta_arc));
    y_arc_outer = r_band_outer * sin(fliplr(theta_arc));

    % Concatenate inner and outer arc points
    x_patch = [x_arc_inner, x_arc_outer];
    y_patch = [y_arc_inner, y_arc_outer];

    % Plot each arc with corresponding color
    patch(x_patch, y_patch, rcolor(k,:), 'EdgeColor', 'none');

end
