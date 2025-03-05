RedCell=RedSignal>=RedTh;
CellColors=repmat([0 1 0],length(cellPlane),1);
CellColors(RedCell,:)=repmat([1 0 0],sum(RedCell),1);


PlotParam.RowPlot=1;
PlotParam.RowColNum=1;
PlotParam.RowColID=1;
PlotParam.EdgeParam=[0.06 0.1 0.1 0.06 0.06 0.06];
PlotParam.CellCenterWith=1;
PlotParam.CellBoundaryWidth=1;
PlotParam.PlotCenter=0;

figure;      
ImgClim=[0 0.9]
H=MultiPlanes2DShow(RegGreen, cellBoundary, Pos3D, [], PlaneZ, CellColors, ImgClim,PlotParam);
b=colorbar;
set(b,'position',[0.95 0.3 0.01 0.5],'ticks',ImgClim);b.Label.String='Amplitude RedCh';
papersizePX=[0 0 21 7]*1.5;
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
% saveas(gcf,[SavePath 'Red'],'fig')
% saveas(gcf,[SavePath 'Red.png'],'png')

figure;      
H=MultiPlanes2DShow(permute(RegRed, [1, 2, 3]), cellBoundary, Pos3D, [], PlaneZ, CellColors, ImgClim,PlotParam);
b=colorbar;
set(b,'position',[0.95 0.3 0.01 0.5],'ticks',ImgClim);b.Label.String='Amplitude GreenCh';
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
% saveas(gcf,[SavePath 'Green'],'fig')
% saveas(gcf,[SavePath 'Green.png'],'png')

