function   [SLMResTemp, mergedRegion, contourPixels] = ROIcontour_center(roiHeatAve,ROIparam)

CellSize=ROIparam.CellSize;  %%by um;
threshold_percentage=ROIparam.threshold_percentage;
max_distance=ROIparam.max_distance;  %% 2/3 diameter of a cell by pixel as maximal response region-SLM center distance
min_region_size=ROIparam.min_region_size;
PeakTh=ROIparam.PeakTh;
min_merged_region_size=ROIparam.min_merged_region_size;
contourMethod=ROIparam.contourMethod;
NeighbourHfWidthPixel=ROIparam.NeighbourHfWidthPixel;

   P.xLeft=0.06;        %%%%%%Left Margin
   P.xRight=0.02;       %%%%%%Right Margin
   P.yTop=0.02;         %%%%%%Top Margin
   P.yBottom=0.06;      %%%%%%Bottom Margin
   P.xInt=0.01;         %%%%%%Width-interval between subplots
   P.yInt=0.005;         %%%%%%Height-interval between subplots


NeighbourRange=[-1 1]*NeighbourHfWidthPixel;


[SLMResTemp, mergedRegion, contourPixels] = check_high_value_center(roiHeatAve, threshold_percentage, max_distance, min_region_size, PeakTh, min_merged_region_size, contourMethod);
