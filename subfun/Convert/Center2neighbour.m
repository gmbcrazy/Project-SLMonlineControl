function roiNeighbour=Center2neighbour(MeanFieldMap,MedCenter,PixFromMedCenter)

%% ROI2neighbour - Extracts the neighborhood of specified ROIs in an image.
% [roiNeighbour, ROIboundary] = ROI2neighbour(MeanFieldMap, cellIDMap, cellIDNeed, PixFromMedCenter)
% 
% Inputs:
%   - MeanFieldMap: Mean field matrix of imaging data.
%   - cellIDMap: Matrix where ROIs are marked with different cell IDs.
%   - cellIDNeed: Cell IDs for which to extract the ROI neighborhood.
%   - PixFromMedCenter: Number of pixels from the median center of each ROI to get the ROI neighborhood.
%
% Outputs:
%   - roiNeighbour: 3D matrix containing the neighborhood around specified ROIs.
%   - ROIboundary: Cell array where each cell contains the boundary coordinates of a specified ROI.



[m,n]=size(MeanFieldMap);
cellIDMapExt=zeros(m+2*PixFromMedCenter,n+2*PixFromMedCenter)+nan;
% cellIDMapExt([1:m]+PixFromMedCenter,[1:n]+PixFromMedCenter)=cellIDMap;


% Extend the original imaging by 2*PixFromMedCenter of pixels, where the original imaging locates in the center part of the extended ones.
MeanFieldMapExt=zeros(m+2*PixFromMedCenter,n+2*PixFromMedCenter)+nan;
MeanFieldMapExt([1:m]+PixFromMedCenter,[1:n]+PixFromMedCenter)=MeanFieldMap;

ROIboundary={};

roiNeighbour=zeros(PixFromMedCenter*2+1,PixFromMedCenter*2+1,size(MedCenter,1))+nan;

for iC=1:size(MedCenter,1)
   
    MedCenter(iC,1)=max(min(MedCenter(iC,1),m),1)+PixFromMedCenter;  %+PixFromMedCenter: noted that the imaging is extented PixFromMedCenter pixels at both head and tails
    MedCenter(iC,2)=max(min(MedCenter(iC,2),n),1)+PixFromMedCenter;  %+PixFromMedCenter: noted that the imaging is extented PixFromMedCenter pixels at both head and tails

end

% Loop through specified cell IDs to extract ROI neighborhood
for iC=1:size(MedCenter,1)
    % [ypix,xpix]=find(cellIDMap==cellIDNeed(iC));    

     % Calculate the median center of the current ROI
    % MedCenter(iC,1)=max(min(MedCenter(iC,1),m),1)+PixFromMedCenter;  %+PixFromMedCenter: noted that the imaging is extented PixFromMedCenter pixels at both head and tails
    % MedCenter(iC,2)=max(min(MedCenter(iC,2),n),1)+PixFromMedCenter;  %+PixFromMedCenter: noted that the imaging is extented PixFromMedCenter pixels at both head and tails

     % Define the range for extracting the ROI neighborhood
    yrange=MedCenter(iC,1)-PixFromMedCenter:MedCenter(iC,1)+PixFromMedCenter;
    xrange=MedCenter(iC,2)-PixFromMedCenter:MedCenter(iC,2)+PixFromMedCenter;

    % Extract the ROI neighborhood and adjust its intensity
    % roiNeighbour(:,:,iC)=imadjust(MeanFieldMapExt(yrange,xrange));  
    roiNeighbour(:,:,iC)=MeanFieldMapExt(yrange,xrange);  

    % Extract the ROI neighborhood IDs
    % temp=cellIDMapExt(yrange,xrange);
    % temp(temp~=cellIDNeed(iC))=0;
    % roiNeighbourID(:,:,iC)=temp;

    
end

