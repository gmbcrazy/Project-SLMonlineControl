function PlotImagePlane(MeanImgPlane)
% This function visualizes targets on multiple image planes, with targets grouped by GroupID.

% Determine the dimensions of the input image
Dim = size(MeanImgPlane);
if length(Dim)==2
   % Single plane image
   PlaneN=1;
else
   % Multi-plane image
   PlaneN=Dim(3);  % Number of planes
end

% Loop through each image plane
for i=1:PlaneN
   % Create subplot for each plane
   subplot(PlaneN,1,i)
   % Normalize amplitude of the image slice and display it
   % tempImg=AmpNormalize(squeeze(MeanImgPlane(:,:,i)),[1 99]);
   tempImg=squeeze(MeanImgPlane(:,:,i));
  
   imagesc(tempImg);
   colormap("gray")  % Use grayscale colormap
   caxis([0 1])  % Set color axis scaling
   ylabel(['Plane ' num2str(i)])  % Label y-axis with plane number
   set(gca,'xtick',[],'ytick',[])  % Remove tick marks


   % Label the x-axis of the last subplot with image dimensions
   if i==PlaneN
       xlabel([num2str(size(tempImg,1)) ' x ' num2str(size(tempImg,1)) ' pixels'])
   end
end
