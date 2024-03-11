
function ROIPairPlot(ROIInfo)

   P.xLeft=0;        %%%%%%Left Margin
   P.xRight=0;       %%%%%%Right Margin
   P.yTop=0;         %%%%%%Top Margin
   P.yBottom=0;      %%%%%%Bottom Margin
   P.xInt=0.000;         %%%%%%Width-interval between subplots
   P.yInt=0.000;         %%%%%%Height-interval between subplots

for i=1:length(ROIInfo.FixRoiB)
    subplotLU(2,length(ROIInfo.FixRoiB),1,i,P)
    imagesc(ROIInfo.FixRoiN(:,:,i));hold on;
    for iB=1:length(ROIInfo.FixRoiB{i})
        plot(ROIInfo.FixRoiB{i}{iB}(:,2),ROIInfo.FixRoiB{i}{iB}(:,1),'r','LineWidth',2);
    end
    if i==1
       ylabel('Previous Day')
    end
    set(gca,'xtick',[],'ytick',[])

    % subplotLU(3,10,2,i)
    % imagesc(MovingRoiRegN(:,:,i));hold on;
    % for iB=1:length(MovingRoiRegB{i})
    %     plot(MovingRoiRegB{i}{iB}(:,2),MovingRoiRegB{i}{iB}(:,1),'r','LineWidth',2);
    % end


    subplotLU(2,length(ROIInfo.FixRoiB),2,i,P)
    imagesc(ROIInfo.MovingRoiN(:,:,i));hold on;
    for iB=1:length(ROIInfo.MovingRoiB{i})
        plot(ROIInfo.MovingRoiB{i}{iB}(:,2),ROIInfo.MovingRoiB{i}{iB}(:,1),'r','LineWidth',2);
    end

    text(20,20,showNum(ROIInfo.corrMatch(i),2),'color','r')

    if i==1
       ylabel('Current Day')
    end
    set(gca,'xtick',[],'ytick',[])


end
