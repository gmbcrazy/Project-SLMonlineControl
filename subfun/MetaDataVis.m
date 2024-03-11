function  MetaDataVis(Data,SaveP)

if nargin==1
   SaveP=0;
end
load(Data);
tempFile=dir(Data);

h = figure('units','normalized','outerposition',[0 0 1 1]);

subplot(10,1,[1 5])
climits = [0 10]; % set colorbar limits to be appropriate
iscell = iscell; ce = find(iscell(:,1) == 1); fDeltaFoF1 = fDeltaFoF(:,ce);
imagesc(fDeltaFoF1',climits)
c = colorbar;
pos = get(c,'Position');
c.Position = [pos(1)+.04, pos(2), 0.01, pos(4)]; % move colorbar so graphs are aligned
c.Label.String = 'dF Z-score)';
colormap bone
title('\fontsize{11}dF')
ylabel('Cell #')
xticks(3600);
xticklabels('2 min')

subplot(10,1,6)
plot(fWhisk1)
hold on
plot(fWhisk2)
xlim([1 size(fWhisk1,1)]);
ylim([0 max(fWhisk1)+0.3]);
title('\fontsize{11}Whisking')
box off
xticks(3600);
xticklabels('\fontsize{10}2 min')



subplot(10,1,7)
if ~isempty(fSpeed)
plot(fSpeed, 'Color', [.85 .33 .1], 'LineWidth', 1)
xlim([1 size(fSpeed,1)]);
ylim([0 max(fSpeed)+0.3]);
end
title('\fontsize{11}Locomotion')
ylabel('\fontsize{10}Speed')
box off
xticks(3600);
xticklabels('\fontsize{10}2 min')
xlim([1 size(fSpeed,1)]);
ylim([0 max(fSpeed)+0.3]);

subplot(10,1,8)
if ~isempty(fStim)
plot(fStim/max(fStim), 'Color', [1 .07 .65], 'LineWidth', 1)
hold on
plot(fSound/max(fSound), 'Color', [0.3 .075 .93], 'LineWidth', 1)
plot(fWhisk1Orig/max(fWhisk1Orig),'Color',[0.1 0.7 0.1],'LineWidth',1)

xlim([1 size(fStim,1)]);
ylim([0 1+0.3]);
title('\fontsize{11}Whisker stim (pink) OR Stimulator Sound (purple) OR RawWhisking + Artifact (green)')
box off
xticks(3600);
xticklabels('\fontsize{10}2 min')
end

subplot(10,1,[9 10])
% plot(fSpeed, 'Color', [.85 .33 .1], 'LineWidth', 1)
% xlim([1 size(fSpeed,1)]);
% ylim([0 max(fSpeed)+0.3]);
% title('\fontsize{11}Locomotion')
% ylabel('\fontsize{10}Speed')
% box off
% xticks(3600);
% xticklabels('\fontsize{10}2 min')
% a=nanzscore(fBehav());
% imagesc(a')
if ~isempty(fBehav)
a=fBehav(:,1:end-1)';
a=nanzscore(a')';
% plot(a')
imagesc(a);
BehavNum=size(fBehav,2);
set(gca,'ylim',[0 BehavNum-0.5],'ytick',[1:BehavNum-1],'yticklabel',BehavStruc(1).BehLabel(1:BehavNum-1));
b = colorbar;
pos = get(c,'Position');
b.Position = [pos(1), pos(2)-0.3, 0.01, 0.15]; % move colorbar so graphs are aligned
% c.Position = [pos(1)+.04, pos(2), 0.01, pos(4)]; % move colorbar so graphs are aligned

b.Label.String = 'BodyMove Z-score)';
colormap jet
title('\fontsize{11}Body Movement')
box off
xticks(3600);
xticklabels('\fontsize{10}2 min')

end


% subplot(10,1,10)
% plot(fLED/max(fLED), 'Color', [1 .07 .65], 'LineWidth', 1)
% hold on
% plot(fLEDCL, 'Color', [0.3 .075 .93], 'LineWidth', 1)
% xlim([1 size(fLED,1)]);
% ylim([0 1+0.3]);
% title('\fontsize{11}LED open-Loop (pink) OR LED closed-loop (purple)')
% box off
% xticks(3600);
% xticklabels('\fontsize{10}2 min')

% [~, figName, ~] = fileparts([currExpVars.folder '\' currExpVars.name]);
% saveas(h, [currExpVars.folder '\' figName '_Raster'])
tempFile.name;
if SaveP==0
   return
elseif SaveP==1
   SaveFile=[tempFile.folder '\' 'tempFile.name' '_Raster'];
   saveas(h, SaveFile)

elseif isstr(SaveP)
   SaveFile=[SaveP '\' 'tempFile.name' '_Raster'];
   saveas(h, SaveFile)

else
   return
end
