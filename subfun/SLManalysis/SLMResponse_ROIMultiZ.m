function [SLMRes,sampleN]=SLMResponse_ROIMultiZ(ROIMap,SLMTrialInfo,ROIparam,minTrialN,SumDataFolder,FileIDrange)

PointAll=ROIparam.PointAll;
PointsTest=ROIparam.PointsTest;
LaserPower=ROIparam.LaserPower;
if isempty(FileIDrange)
   FileIDrange=[min(SLMTrialInfo(:,1));max(SLMTrialInfo(:,1))];
end
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

% Check exsiting ROI png to avoid overwrite
ExistingPngList=dir([SumDataFolder 'CheckROI*.png']);
FileID=length(ExistingPngList)+1;
OutPng=[SumDataFolder 'CheckROI' num2str(FileID) '.png'];



% CellSize=20;  %%by um;
% threshold_percentage=0.3;
% thNum=15;
% max_distance=ceil(CellSize*1/4/ROIparam.umPerPixel);  %% 1/4 diameter of a cell by pixel as maximal response region-SLM center distance
% min_region_size=10;
% PeakTh=200;
% min_merged_region_size=20;
% contourMethod='perim';

SLMRes=zeros(length(PointAll),length(ROIparam.LaserPower));
sampleN=SLMRes;

figure;
for iPoint=1:length(PointsTest)
    for iLaser=1:length(LaserPower)
        Point=PointsTest(iPoint);
        I1=find((SLMTrialInfo.Point==Point&ismember(SLMTrialInfo.Laser,LaserPower(iLaser)))==1);
        % limit range FileID in binFile 
        I2=[];
        for jFile=1:size(FileIDrange,2)
        I2=union(I2,find(SLMTrialInfo.FileID<=FileIDrange(2,jFile)&SLMTrialInfo.FileID>=FileIDrange(1,jFile)));
        end
        I1=intersect(I1,I2);
        sampleN(Point,iLaser)=length(I1);
        if length(I1)>1
           roiHeatAve=squeeze(mean(ROIMap(:,:,I1,:),3));  
        elseif length(I1)==1
           roiHeatAve=squeeze(ROIMap(:,:,I1,:));  
        else
           continue
        end
        [SLMResTemp, mergedRegion, contourPixels] = check_high_value_center(roiHeatAve, threshold_percentage, max_distance, min_region_size, PeakTh, min_merged_region_size, contourMethod);
        SLMRes(Point,iLaser)=SLMResTemp;

        subplotLU(length(PointAll),length(LaserPower),Point,iLaser,P);


        imagesc(roiHeatAve);hold on;
        if ~isempty(mergedRegion)&&SLMResTemp
            plot(contourPixels(:,2),contourPixels(:,1),'g.');
        end
        set(gca,'clim',ROIparam.Clim,'xtick',[],'ytick',[]);
        colormap(ROIparam.Colormap);
        % text(NeighbourHfWidthPixel,0,['P' num2str(Point) ', ' num2str(length(I1))],'color',[0 0 0]) 

        if length(I1)>=minTrialN
           if SLMRes(Point,iLaser)
            text(NeighbourHfWidthPixel,0,['P' num2str(Point) ', ' num2str(length(I1))],'color',[0 1 0]); 
           else
            text(NeighbourHfWidthPixel,0,['P' num2str(Point) ', ' num2str(length(I1))],'color',[0 0 0]); 
           end
        else
            text(NeighbourHfWidthPixel,0,['P' num2str(Point) ', ' num2str(length(I1))],'color',[0 0 0]);
        end




        % % if length(I1)>=minTrialN         
        % %    subplotLU(length(PointsAll),length(LaserPower),Point,iLaser,P);
        % %    % PlaneI=find(abs(ROIparam.TotalSLMPos3D(Point,3)-ROIparam.PlaneZ)<1);
        % %    roiHeatAve=squeeze(mean(ROIMap(:,:,I1,:),3));  
        % %    % heatMap=squeeze(PSTHtemp(:,:,PlaneI))';
        % %    % roiHeat=Center2neighbour(heatMap,ROIparam.TotalSLMPos3D(Point,[2,1]),NeighbourHfWidthPixel);
        % %    [SLMRes(Point,iLaser), mergedRegion, contourPixels] = check_high_value_center(roiHeatAve, threshold_percentage, max_distance, min_region_size, PeakTh, min_merged_region_size, contourMethod);
        % % 
        % %    sampleN(Point,iLaser)=length(I1);
        % %    imagesc(roiHeatAve);hold on;
        % %    set(gca,'clim',ROIparam.Clim,'xtick',[],'ytick',[])
        % %    colormap(ROIparam.Colormap)
        % % % for iField=1:length(regions)
        % %    if ~isempty(mergedRegion)
        % %     plot(contourPixels(:,2),contourPixels(:,1),'g.');
        % %    end
        % %  text(NeighbourHfWidthPixel,0,['P' num2str(Point) ', ' num2str(length(I1))],'color',[0 1 0]) 
        % % end


    end
end

papersizePX=[0 0 length(LaserPower)*4 length(PointsTest)*4];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
saveas(gcf,OutPng,'png');
% close all
