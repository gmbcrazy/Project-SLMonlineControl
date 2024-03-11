function [imgStack] = readTifsAna2019(dirname, h, w, specialInstruc)

% specialInstruc is optional parameter in char form. Example is 'Ch2'

if nargin == 3
    imagefiles = dir([dirname '\*.tif']);
    nfiles = length(imagefiles);
else
    imagefiles = dir([dirname '\*' specialInstruc '*.tif']);
    nfiles = length(imagefiles);
end
% Assume all files are 1-page tifs
images = uint16(zeros(h, w, nfiles)); % change this to fit img resolution

for i = 1:nfiles
    currentfilename = [dirname '\' imagefiles(i).name];
    image = imread(currentfilename);
    images(:,:,i) = image;
end 
imgStack = images;

disp(['Finished loading ' num2str(size(imgStack,3)) ' images.']);
