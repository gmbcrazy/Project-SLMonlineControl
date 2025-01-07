function [SLMRes,sampleN,DeltaRatio,p]=SLMResponse_PSTHTraceMultiZ(PSTHTrace,SLMTrialInfo,PSTHparam,minTrialN,SumDataFolder,FileIDrange)

PointAll=PSTHparam.PointAll;
PointsTest=PSTHparam.PointsTest;
LaserPower=PSTHparam.LaserPower;
if isempty(FileIDrange)
   FileIDrange=[min(SLMTrialInfo(:,1));max(SLMTrialInfo(:,1))];
end
% CellSize=PSTHparam.CellSize;  %%by um;
% threshold_percentage=PSTHparam.threshold_percentage;
% max_distance=PSTHparam.max_distance;  %% 2/3 diameter of a cell by pixel as maximal response region-SLM center distance
% min_region_size=PSTHparam.min_region_size;
% PeakTh=PSTHparam.PeakTh;
% min_merged_region_size=PSTHparam.min_merged_region_size;
% contourMethod=PSTHparam.contourMethod;
% NeighbourHfWidthPixel=PSTHparam.NeighbourHfWidthPixel;

   P.xLeft=0.06;        %%%%%%Left Margin
   P.xRight=0.02;       %%%%%%Right Margin
   P.yTop=0.02;         %%%%%%Top Margin
   P.yBottom=0.06;      %%%%%%Bottom Margin
   P.xInt=0.01;         %%%%%%Width-interval between subplots
   P.yInt=0.005;         %%%%%%Height-interval between subplots



% Check exsiting ROI png to avoid overwrite
ExistingPngList=dir([SumDataFolder 'CheckTrace*.png']);
FileID=length(ExistingPngList)+1;
OutPng=[SumDataFolder 'CheckTrace' num2str(FileID) '.png'];


% PSTHparam.Clim=[-500 500]

% CellSize=20;  %%by um;
% threshold_percentage=0.3;
% thNum=15;
% max_distance=ceil(CellSize*1/4/PSTHparam.umPerPixel);  %% 1/4 diameter of a cell by pixel as maximal response region-SLM center distance
% min_region_size=10;
% PeakTh=200;
% min_merged_region_size=20;
% contourMethod='perim';

SLMRes=zeros(length(PointAll),length(PSTHparam.LaserPower));
DeltaRatio=SLMRes+nan;
p=SLMRes+nan;

sampleN=SLMRes;
X=1:PSTHparam.PreSLMCal+PSTHparam.PostSLMCal;
BaseI=1:PSTHparam.PreSLMCal;
ResI=PSTHparam.PreSLMCal+1:PSTHparam.PreSLMCal+PSTHparam.PostSLMCal;




figure;
for iPoint=1:length(PointsTest)
    for iLaser=1:length(LaserPower)
        Point=PointsTest(iPoint);
        % I1=find((SLMTrialInfo.Point==Point&ismember(SLMTrialInfo.Laser,LaserPower(iLaser)))==1);
        I1=find((SLMTrialInfo.Point==Point&abs(SLMTrialInfo.Laser-LaserPower(iLaser))<0.01)==1);

        % limit range FileID in binFile 
        I2=[];
        for jFile=1:size(FileIDrange,2)
        I2=union(I2,find(SLMTrialInfo.FileID<=FileIDrange(2,jFile)&SLMTrialInfo.FileID>=FileIDrange(1,jFile)));
        end
        I1=intersect(I1,I2);
        sampleN(Point,iLaser)=length(I1);
        if length(I1)>1
           traceAve=mean(PSTHTrace(I1,:),1);
           traceSte=ste(PSTHTrace(I1,:));
        elseif length(I1)==1
           traceAve=PSTHTrace(I1,:); 
           traceSte=zeros(size(traceAve));
        else
           continue
        end
        subplotLU(length(PointAll),length(LaserPower),Point,iLaser,P);
        % plot(PSTHTrace(I1,:)','Color',[0.1 0.1 0]);
        error_area(X,traceAve,traceSte,[0.1 0.1 0.1],0.4);
        text(5,200,['P' num2str(Point) ', ' num2str(length(I1))],'color',[0 0 0]);
        set(gca,'xlim',[0 max(X)],'ylim',PSTHparam.YLim,'ytick',[0 PSTHparam.YLim(end)]);
        plot([PSTHparam.PreSLMCal PSTHparam.PreSLMCal]+0.5,PSTHparam.YLim,':','Color',[0 0 0]);
        
        if length(I1)>=minTrialN         
           [DeltaRatio(Point,iLaser), p(Point,iLaser)] = PSTHtrace_PostPreTest(PSTHTrace(I1,:), BaseI, ResI,PSTHparam.TestMethod);
           if DeltaRatio(Point,iLaser)>0&&p(Point,iLaser)<=PSTHparam.pTh
              SLMRes(Point,iLaser)=1;
              text(5,400,['DeltaR: ' sprintf('%.2f', DeltaRatio(Point,iLaser)) ', Pvalue: ' sprintf('%.3f', p(Point,iLaser))],'color',[0 1 0]);
           end
        end
        % if length(I1)>=minTrialN         
        %    subplotLU(length(PointsTest),length(LaserPower),iPoint,iLaser,P);
        %    % PlaneI=find(abs(PSTHparam.TotalSLMPos3D(Point,3)-PSTHparam.PlaneZ)<1);
        %    traceAve=mean(PSTHTrace(I1,:),1);  
        %    % heatMap=squeeze(PSTHtemp(:,:,PlaneI))';
        %    % roiHeat=Center2neighbour(heatMap,PSTHparam.TotalSLMPos3D(Point,[2,1]),NeighbourHfWidthPixel);
        %    % [SLMRes(Point,iLaser), mergedRegion, contourPixels] = check_high_value_center(roiHeatAve, threshold_percentage, max_distance, min_region_size, PeakTh, min_merged_region_size, contourMethod);
        % 
        % 
        %    subplotLU(length(PointsAll),length(LaserPower),Point,iLaser,P);
        %    % a=error_area(1:length(traceAve),traceAve,ste(PSTHTrace),[0 1 0],0.5);
        % 
        %    % set(gca,'clim',PSTHparam.Clim,'xtick',[],'ytick',[])
        %    % colormap(PSTHparam.Colormap)
        % % for iField=1:length(regions)
        % % if ~isempty(mergedRegion)
        % %     plot(contourPixels(:,2),contourPixels(:,1),'g.');
        % % end
        %    text(10,50,['P' num2str(Point) ', ' num2str(length(I1))],'color',[0 1 0]) 
        % end


    end
end

papersizePX=[0 0 length(LaserPower)*4 length(PointsTest)*4];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
saveas(gcf,OutPng,'png');
% close all
