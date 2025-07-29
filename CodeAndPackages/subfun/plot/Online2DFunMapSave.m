
function Online2DFunMapSave(CaData,cellBoundary, Pos3D, PlaneZ, CellSpeedColors, ImgClim, PlotParam, SavePathAllPoint, BehLabel)

figure;      
H=MultiPlanes2DShow(permute(CaData.PlaneMeanImg, [2, 1, 3]), cellBoundary, Pos3D, [], PlaneZ, CellSpeedColors, ImgClim,PlotParam);
subplot('position',[0.95 0.3 0.01 0.3])
MultibarPlot(1:64,gray(64),0,0,1,1);
set(gca,'xlim',[0 64],'xtick',[0 32 64],'xticklabel',sort([ImgClim 0]),'ytick',[]);
xlabel(BehLabel)
camroll(90);
papersizePX=[0 0 10*length(PlaneZ) 10];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
saveas(gcf,[SavePathAllPoint BehLabel ],'fig')
saveas(gcf,[SavePathAllPoint BehLabel '.png'],'png')
print(gcf,[SavePathAllPoint BehLabel '.svg'], '-dsvg', '-painters');

close all
