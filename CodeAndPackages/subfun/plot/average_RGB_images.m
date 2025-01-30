% Script to take the gray scale or RGB images in the folder and average them together to produce a mean image.
% If the images are of different sizes, it resizes them to the size of the first image.
% If the images are of different color bits, it makes them the color (gray scale or RGB) of the first image.
clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
clear;  % Erase all existing variables. Or clearvars if you want.
workspace;  % Make sure the workspace panel is showing.
format long g;
format compact;
fontSize = 20;

%--------------------------------------------------------------------------------
% Get folder name where the images are that we want to average.
% Option #1:
folder = fileparts(which('cameraman.tif')); % Determine where MATLAB demo images folder is.
% Option #2:
folder = 'D:\My Pictures\Illusions'; % Specify some particular folder.
% Option #3:
folder = uigetdir(pwd, 'Select folder');
% folder will be 0 (a double) if they click cancel.
% folder will be the path (a string) if they clicked OK.
if folder == 0
	% Clicked cancel.  Exit program.
	return;
end
% Comment out whichever folder selection options above that you don't want to use.

% Make sure the folder actually exists.
if ~isdir(folder)
	errorMessage = sprintf('Error: The following folder does not exist:\n%s', folder);
	uiwait(warndlg(errorMessage));
	return;	
else
	fprintf('Averaging images in the following folder:\n          %s', folder);
end

%--------------------------------------------------------------------------------
% Get a list of all TIF, PNG, BMP, and JPG files.
imageFiles = [dir(fullfile(folder,'*.TIF')); dir(fullfile(folder,'*.PNG')); dir(fullfile(folder,'*.BMP')); dir(fullfile(folder,'*.jpg'))];

%--------------------------------------------------------------------------------
% Now do the averaging in a loop
numberOfImages = length(imageFiles);
theyreColorImages = false;
for k = 1 : numberOfImages
	fullFileName = fullfile(folder, imageFiles(k).name);
	fprintf('About to read %s\n', fullFileName);
	thisImage=imread(fullFileName);
	[thisRows, thisColumns, thisNumberOfColorChannels] = size(thisImage);
	if k == 1
		% Save the first image.
		sumImage = double(thisImage);
		% Save its dimensions so we can match later images' sizes to this first one.
		rows1 = thisRows;
		columns1 = thisColumns;
		numberOfColorChannels1 = thisNumberOfColorChannels;
		theyreColorImages = numberOfColorChannels1 >= 3;

		% Add in the "if" block below if you want to force them to be color, even if the first image is gray scale.
% 		if numberOfColorChannels1 == 1
% 			% Convert to color
% 			thisImage = cat(3, thisImage, thisImage, thisImage);
% 			theyreColorImages = false;
% 		end
	else
		% It's the second, or later, image.
		if rows1 ~= thisRows || columns1 ~= thisColumns
			% It's not the same size, so resize it to the size of the first image.
			thisImage = imresize(thisImage, [rows1, columns1]);
		end
		% Make sure the colors match - either all color or all gray scale, according to the first one.
		if thisNumberOfColorChannels == 3 && numberOfColorChannels1 == 1
			% We have color.  Need to change it to grayscale to match the first one.
			thisImage = rgb2gray(thisImage);
			theyreColorImages = false;
		elseif thisNumberOfColorChannels == 1 && numberOfColorChannels1 == 3
			% We have grayscale.  Need to change it to RGB to match the first one.
			thisImage = cat(3, thisImage, thisImage, thisImage);
			theyreColorImages = true;
		end
		% Now do the summation.
		sumImage = sumImage + double(thisImage); % Be sure to cast to double to prevent clipping.	[rows, columns, numberOfColorBands]=size(thisImage);
		
		% It can't display an RGB image if it's floating point and more than 255.
		% So divide it by the number of images to get it into the 0-255 range.
		if theyreColorImages
			displayedImage = uint8(sumImage / k);
		else
			displayedImage = sumImage;
		end
		imshow(displayedImage, []);
		drawnow;
	end
end

%--------------------------------------------------------------------------------
% Compute and display the final image:
sumImage = uint8(sumImage / numberOfImages);
cla;
imshow(sumImage, []);
caption = sprintf('Average of %d Images', numberOfImages);
title(caption, 'FontSize', fontSize);

%--------------------------------------------------------------------------------
% Set up figure properties:
% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% Get rid of tool bar and pulldown menus that are along top of figure.
set(gcf, 'Toolbar', 'none', 'Menu', 'none');
% Give a name to the title bar.
set(gcf, 'Name', 'Demo by ImageAnalyst', 'NumberTitle', 'Off') 

