% TestHologram
% Lloyd Russell 2016
% Generate a random SLM phase mask to test import and use of HOTlab
% library. Implements a weighted GS algorithm using Martin Perssons HOTlab
% (https://github.com/MartinPersson/HOTlab)


% Load HOTlab
if ~libisloaded('GenerateHologramCUDA')
    [notfound, warnings] = loadlibrary('GenerateHologramCUDA', 'GenerateHologramCUDA.h');
    libfunctions('GenerateHologramCUDA')
end


% Start CUDA
SLMsize = 512;
deviceId = 0;
SLMphasemask = zeros(SLMsize,SLMsize);  % the starting phase mask (seed)
err1 = calllib('GenerateHologramCUDA','startCUDA', SLMphasemask, deviceId);


% Generate Hologram
N_spots = 5;
radial_dist_from_zo = 25;
x_spots = randi(2*radial_dist_from_zo,N_spots,1) - radial_dist_from_zo;
y_spots = randi(2*radial_dist_from_zo,N_spots,1) - radial_dist_from_zo;
z_spots = zeros(N_spots);
I_spots = ones(N_spots);
N_iterations = 30;
method = 2;  % 0: Complex addition of "Lenses and Prisms", no optimization (3D)
             % 1: Weighted Gerchberg-Saxton algorithm using Fresnel propagation (3D)
             % 2: Weighted Gerchberg-Saxton algorithm using fast fourier transforms (2D)
tic
[err2,~,SLMphasemask,~,~,~,~,~] = calllib('GenerateHologramCUDA','GenerateHologram',...
    [], SLMphasemask, x_spots, y_spots, z_spots, I_spots, N_spots, N_iterations, [], method);
toc


% Stop CUDA
err3 = calllib('GenerateHologramCUDA','stopCUDA');
unloadlibrary('GenerateHologramCUDA');  % tidy up, unload dll


% Plot
figure('Position', [100 100 600 300])
subplot(1,2,1)
plot(x_spots,y_spots, 'k.')
axis square
xlim([-SLMsize/2 SLMsize/2])
ylim([-SLMsize/2 SLMsize/2])
title('SLM spot targets')

subplot(1,2,2)
imshow(SLMphasemask)
title('SLM phase mask')