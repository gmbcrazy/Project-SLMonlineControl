function theta=ClockWiseNode(N)

% Get the number of nodes
N = N;

% Define the circle radius
% radius = 10; % Adjust as needed

% Compute angles (starting from 12 o'clock position)
angles = linspace(pi/2, -3*pi/2, N+1);  % Starts at 12 o'clock and goes clockwise

% Compute new X, Y coordinates
% X_new = radius * cos(angles(1:N));
% Y_new = radius * sin(angles(1:N));

theta=angles(1:N);
