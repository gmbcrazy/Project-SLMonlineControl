
function theta=ClockWiseGraph1(G,p,radius)

% Get the number of nodes
N = size(G.Nodes, 1)-1;

% Define the circle radius
% radius = 10; % Adjust as needed

% Compute angles (starting from 12 o'clock position)
angles = linspace(pi/2, -3*pi/2, N+1);  % Starts at 12 o'clock and goes clockwise

% Compute new X, Y coordinates
X_new = radius * cos(angles(1:N));
Y_new = radius * sin(angles(1:N));

theta=angles(1:N);

% Update graph layout
p.XData = [X_new 0];
p.YData = [Y_new 0];

end
