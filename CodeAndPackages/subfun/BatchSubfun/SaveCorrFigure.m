function SaveCorrFigure(NData, rStim, BehavTseries, InitialInd, Nlabel, ResultFolder,BehLabel)
% Plot and save sorted stimulus correlation heatmap/trace for given data type

[rStimSorted, rankStimTemp] = sort(rStim, 'descend');
figure;
subplot('position',[0.1 0.1 0.7 0.6])
imagesc(AmpNormalizeRow(NData(rankStimTemp, InitialInd), [1 99]))
set(gca, 'xlim', [0 length(InitialInd)], 'ylim', [0 size(NData,1)], 'clim', [0.1 0.9])
xlabel('Time (frames)')
ylabel('Cells')
colormap(gray)
b=colorbar;
bar1=colorbar;
bar1.Location = 'northoutside';
bar1.Position = [0.4 0.92 0.2 0.01];
bar1.Label.String = Nlabel;

subplot('position',[0.1 0.75 0.7 0.1])
plot(BehavTseries)
set(gca,'xlim',[0 length(InitialInd)],'xticklabel',[],'Box','off');
ylabel(BehLabel)

subplot('position',[0.85 0.1 0.1 0.6])
plot(rStimSorted,'color',[0.01 0.01 0.1]);
view(90, 90);
set(gca,'xticklabel',[],'Box','off','xlim',[0 size(NData,1)]);
ylabel(['r' BehLabel ' Corr.'])
papersizePX=[0 0 16 16];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
saveas(gcf,[ResultFolder 'Spon' BehLabel Nlabel 'Corr.png']);
close(gcf);
end
