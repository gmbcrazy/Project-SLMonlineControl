function PlotSLMTargets(MeanImgPlane,TargetXY_Zplane,GroupID,GroupColor,TargetRadius)
% This function visualizes targets on multiple image planes, with targets grouped by GroupID.

% Determine the dimensions of the input image
Dim = size(MeanImgPlane);
if length(Dim)==2
   % Single plane image
   PlaneN=1;
   CellPlaneID = ones(length(GroupID),1);  % Assign all targets to the first plane
else
   % Multi-plane image
   PlaneN=Dim(3);  % Number of planes
   CellPlaneID=TargetXY_Zplane(:,3);  % Z-plane identifier for each target
end

% Extract XY coordinates and swap X and Y to align with image coordinates
XY=TargetXY_Zplane(:,[2,1]);

% Identify unique groups
GroupCandi=unique(GroupID);

% Loop through each image plane
for i=1:PlaneN
   % Create subplot for each plane
   subplot(PlaneN,1,i)
   % Normalize amplitude of the image slice and display it
   tempImg=AmpNormalize(squeeze(MeanImgPlane(:,:,i)),[1 99]);
   imagesc(tempImg);
   colormap("gray")  % Use grayscale colormap
   caxis([0 1])  % Set color axis scaling
   ylabel(['Plane ' num2str(i)])  % Label y-axis with plane number
   set(gca,'xtick',[],'ytick',[])  % Remove tick marks

   % Loop through each group and plot centers
   for iG=1:length(GroupCandi)
       % Find indices of cells belonging to the current group in the current plane
       CellI=find(GroupID==GroupCandi(iG)&TargetXY_Zplane(:,3)==i);
       hold on;
       % Plot the cell centers
       plotCellCenter(XY(CellI,:),TargetRadius,GroupColor(iG,:),2)
   end

   % Label the x-axis of the last subplot with image dimensions
   if i==PlaneN
       xlabel([num2str(size(tempImg,1)) ' x ' num2str(size(tempImg,1)) ' pixels'])
   end
end
