
function ROIPlot(ROImap,ROIboundary)


    imagesc(ROImap);hold on;
    for iB=1:length(ROIboundary)
        plot(ROIboundary{iB}(:,2),ROIboundary{iB}(:,1),'r','LineWidth',2);
    end

    set(gca,'xtick',[],'ytick',[])




