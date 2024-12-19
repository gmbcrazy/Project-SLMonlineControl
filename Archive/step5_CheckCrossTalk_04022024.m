clear all
ProcessFolder='F:\LuSLMOnlineTest\04022024\DataDelete1stFrame\';

load('C:\Users\zhangl33\Projects\GenMatCode\Plotfun\Color\colorMapPN3.mat');

ResultFolder=[ProcessFolder 'Result\'];
mkdir(ResultFolder)

TiffFolder=dir([ProcessFolder 'TSeries*'])
temp=struct2table(TiffFolder);
TiffFolder=TiffFolder(temp.isdir==1);
for i=1:size(TiffFolder,1)
    Session(i)=str2num(TiffFolder(i,1).name(end-2:end));
end

FrameNum=[4500;ones(length(TiffFolder)-1,1)*39];    %%For deleting the 1st frame data, with missing data; Not the 001 file of 4500 repetitions.


clear SessInfoNeed
SessInfoNeed.Repetitions=FrameNum;
LastFrame=cumsum(FrameNum);
FirstFrame=[1;LastFrame(1:end-1)+1];

Fall=load('F:\LuSLMOnlineTest\04022024\DataDelete1stFrame\suite2p\combined\Fall.mat');
cellInfo=Suite2pCellInfo(Fall);

PreImgN=8;
PostImgN=35;
clear IndStart IndEnd
for i=1:length(FirstFrame)
    IndStart(i,1)=FirstFrame(i)-PreImgN;
    IndEnd(i,1)=FirstFrame(i)+PostImgN-1;
end
SessInfoNeed.PreFrame=IndStart;
SessInfoNeed.PostFrame=IndEnd;

Invalid=SessInfoNeed.PreFrame<0;
SessInfoNeed.PreFrame(Invalid)=[];
SessInfoNeed.PostFrame(Invalid)=[];
SessInfoNeed.Repetitions(Invalid)=[];



% LaserG=unique(SessInfoNeed.LaserPower)
RepG=unique(SessInfoNeed.Repetitions)
% PmtG=unique(SessInfoNeed.PMTLevel)

deltaFoF=F2deltaFoF(Fall.F,Fall.Fneu,Fall.ops.fs);


NData={deltaFoF', Fall.spks};
Nlabel={'DeltaF', 'Spks'}

cellInfo=Suite2pCellInfo(Fall);

planeG=unique(cellInfo.iplane);
for iplane=1:length(planeG)
    PlaneC(iplane)=max(find(cellInfo.iplane==planeG(iplane)));
end


iscell=find(Fall.iscell(:,1)==1);


ClimScale=[-5 5;-10 10]
   P.xLeft=0.2;        %%%%%%Left Margin
   P.xRight=0.2;       %%%%%%Right Margin
   P.yTop=0.02;         %%%%%%Top Margin
   P.yBottom=0.1;      %%%%%%Bottom Margin
   P.xInt=0.01;         %%%%%%Width-interval between subplots
   P.yInt=0.03;         %%%%%%Height-interval between subplots

for iRep=1:1

    I1=find(SessInfoNeed.Repetitions==RepG(iRep));
            figure;

    for iData=1:length(NData)
        % TempData=zscore(double(NData{iData}(iscell,:)),0,2);
        TempData=double(NData{iData}(iscell,:));
        % for iCell=1:size(cellInfo,1)
            % for iLaser=1:length(LaserG)
               
                % I2=find(SessInfoNeed.LaserPower==LaserG(iLaser));
                I3=I1;
                if ~isempty(I3)
                   tempPSTH=[];

                   for iSess = 1:length(I3)
                       TempI=SessInfoNeed.PreFrame(I3(iSess)):SessInfoNeed.PostFrame(I3(iSess));
                       tempPSTH(:,:,iSess)=TempData(:,TempI);
                   end
                   tempPSTH=squeeze(mean(tempPSTH,3));
                   BaseLine=repmat(mean(tempPSTH(:,1:PreImgN),2),1,size(tempPSTH,2));
                   tempPSTH=(tempPSTH-BaseLine)./BaseLine;


                   subplotLU(length(NData),1,iData,1,P);
                   imagesc(tempPSTH);hold on;
                   if iData==2
                   set(gca,'xlim',[0 PreImgN+0.5+PostImgN] ,'xtick',[0 PreImgN+0.5 PreImgN+0.5+PostImgN],'xticklabel',{['-' num2str(PreImgN)] '0' num2str(PostImgN)})
                   % xlabel(['PV-Power ' num2str(LaserG(iLaser))])
                   else
                   set(gca,'xlim',[0 PreImgN+0.5+PostImgN] ,'xtick',[0 PreImgN+0.5 PreImgN+0.5+PostImgN],'xticklabel',[])
                       
                   end
                   set(gca,'tickdir','out','clim',ClimScale(iData,:))
                   colormap(ColorPN3)

                      ylabel('CellsID');
                      set(gca,'ytick',PlaneC+0.5,'yticklabel',{'Plane1' 'Plane2' 'Plane3'})

       

                   for iplane=1:length(PlaneC)
                       plot([0 PreImgN+0.5+PostImgN],zeros(1,2)+PlaneC(iplane)+0.5,'g-')
                   end
                   plot([PreImgN PreImgN]+0.5,[0 length(iscell)],'k-')

                end

 
         b = colorbar;
         set(b,'position',[0.83 0.7-(iData-1)*0.5 0.02 0.2]);
         ylabel(b, Nlabel{iData});

        % end
    end
     papersizePX=[0 0 9 length(NData)*5];
     set(gcf, 'PaperUnits', 'centimeters');
     set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
     saveas(gcf,[ResultFolder 'AllNeuroSessReptition' num2str(RepG(iRep))],'png');
    close all
end


% jet=colormap("jet");
% colorLaser=jet(1:size(jet,1)/length(LaserG):size(jet,1),:);
% close all
colorLaser=[1 0 0];


NData={Fall.F, deltaFoF', Fall.spks};
Nlabel={'F','DeltaF', 'Spks'}

NData={Fall.F, deltaFoF', Fall.spks};
Nlabel={'F','DeltaF', 'Spks'}


for iData=1:length(NData)
figure;
hold on;
    for iCell=1:size(cellInfo,1)

    I1=find(SessInfoNeed.Repetitions==RepG(iRep));
        % TempData=zscore(double(NData{iData}(iscell,:)),0,2);
        TempData=double(NData{iData}(iscell(iCell),:));
        TempData=AmpNormalize(TempData,[0 100]);
        plot(smooth(TempData,3)-iCell*1,'color',[0.2 0.2 0.2]);

        EndPlane=ismember(iCell,PlaneC)
        if EndPlane
           plot([0 length(TempData)],zeros(1,2)-iCell,'g-');
        end

    end
    for iS=1:size(SessInfoNeed.PreFrame,1)
       % [~,Laserj]=ismember(SessInfoNeed.LaserPower(iS),LaserG);
       plot(zeros(1,2)+SessInfoNeed.PreFrame(iS)+PreImgN,[-iCell 0],'color',colorLaser);
       % text(SessInfoNeed.PreFrame(iS)+PreImgN,0,['Laser' num2str(LaserG(Laserj))],'HorizontalAlignment','center');
    end
    set(gca,'ytick',[-length(iscell):1:-1]+0.5,'yticklabel',abs([-length(iscell):1:-1]))
    ylabel('Cell IDs')
     % colormap(colorLaser)    
     % b = colorbar;
     % set(b,'position',[0.9 0.2 0.03 0.5]);
     % ylabel(b, 'PV power');
     % set(b,'ytick',[1:length(LaserG)]/length(LaserG),'yticklabel',LaserG);
     papersizePX=[0 0 20 20];
     set(gcf, 'PaperUnits', 'centimeters');
     set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));

% Assuming TempData and iscell are defined earlier in your script as shown.
xSlider = uicontrol('Style', 'slider', 'Min', 1, 'Max', length(TempData), 'Value', length(TempData)/2, ...
                    'Position', [100 20 300 20], 'Callback', @xSliderCallback);

ySlider = uicontrol('Style', 'slider', 'Min', -length(iscell), 'Max', -1, 'Value', -length(iscell)/2, ...
                    'Position', [20 100 20 300], 'Callback', @ySliderCallback);


     saveas(gcf,[ResultFolder Nlabel{iData} 'AllNeuroSig'],'fig');
     close all
end


ResultFolderCell=[ResultFolder 'Cells\'];
mkdir(ResultFolderCell);

% for iCell=1:size(cellInfo,1)

   P.xLeft=0.1;        %%%%%%Left Margin
   P.xRight=0.1;       %%%%%%Right Margin
   P.yTop=0.02;         %%%%%%Top Margin
   P.yBottom=0.1;      %%%%%%Bottom Margin
   P.xInt=0.01;         %%%%%%Width-interval between subplots
   P.yInt=0.03;         %%%%%%Height-interval between subplots

for iCell=1:size(cellInfo,1)


    figure;
for iRep=1:1
% for iRep=2:2
    I1=find(SessInfoNeed.Repetitions==RepG(iRep));
    for iData=1:length(NData)
        % TempData=zscore(double(NData{iData}(iscell,:)),0,2);
        TempData=double(NData{iData}(iscell(iCell),:));
        TempData=AmpNormalize(TempData,[0 100]);

            % for iLaser=1:length(LaserG)
            % for iLaser=1:8
                % I2=find(SessInfoNeed.LaserPower==LaserG(iLaser));
                I3=I1;
                if ~isempty(I3)
                   tempPSTH=[];
                   for iSess = 2:length(I3)
                       TempI=SessInfoNeed.PreFrame(I3(iSess)):SessInfoNeed.PostFrame(I3(iSess));
                       tempPSTH(:,iSess)=TempData(TempI);
                   end
                   % PSTHLaser(:,iLaser)=squeeze(mean(tempPSTH,2));
                   % error_area(1:size(tempPSTH,1),mean(tempPSTH,2),std(tempPSTH,0,2),colorLaser(iLaser,:),0.5)
                   BaseLine=repmat(mean(tempPSTH(1:PreImgN,:),1),size(tempPSTH,1),1);
                   tempPSTH=tempPSTH-BaseLine;

                   subplotLU(2,length(NData),1,iData,P);hold on
                   error_area(1:size(tempPSTH,1),mean(tempPSTH,2),ste(tempPSTH')',colorLaser,0.5);
                   % if iRep==1&&iLaser==1
                   %    text(15,max(mean(tempPSTH,2)),Nlabel{iData})
                   %    % set(gca,'ylim',[0 0.5])
                   % end

                   % plot(PSTHLaser(:,iLaser),'color',colorLaser(iLaser,:));
                   % set(gca,'xtick',[0 PreImgN+0.5 PreImgN+0.5+PostImgN],'xticklabel',{'Start' num2str(PreImgN) num2str(PreImgN+PostImgN)})
                   % set(gca,'tickdir','out','clim',ClimScale(iData,:))
                   % if iRep==3
                   set(gca,'xlim',[0 PreImgN+0.5+PostImgN] ,'xtick',[0 PreImgN+0.5 PreImgN+0.5+PostImgN],'xticklabel',{['-' num2str(PreImgN)] '0' num2str(PostImgN)})
                   xlabel(Nlabel{iData})
                   set(gca,'ylim',[-0.02 0.12])
                   % else
                   % set(gca,'xlim',[0 PreImgN+0.5+PostImgN] ,'xtick',[0 PreImgN+0.5 PreImgN+0.5+PostImgN],'xticklabel',[])
                   % 
                   % end
                   subplotLU(2,length(NData),2,iData,P);hold on
                   plot(1:size(tempPSTH,1),tempPSTH,'color',[0.5 0.5 0.5]);
                   % if iRep==3
                   set(gca,'xlim',[0 PreImgN+0.5+PostImgN] ,'xtick',[0 PreImgN+0.5 PreImgN+0.5+PostImgN],'xticklabel',{['-' num2str(PreImgN)] '0' num2str(PostImgN)})
                   xlabel(Nlabel{iData})
                   % else
                   % set(gca,'xlim',[0 PreImgN+0.5+PostImgN] ,'xtick',[0 PreImgN+0.5 PreImgN+0.5+PostImgN],'xticklabel',[])
                   % 
                   % end

                   if iData==1
                      ylabel(['Frame #' num2str(RepG(iRep))]);
                   end
            
                   % colormap(ColorPN3)

                end

            % end
     % colormap(colorLaser)    
     % b = colorbar;
     % set(b,'position',[0.92 0.2 0.01 0.5]);
     % ylabel(b, 'PV power');
     % set(b,'ytick',[1:length(LaserG)]/length(LaserG),'yticklabel',LaserG);

    end
end
     papersizePX=[0 0 length(NData)*6 2*5 ];
     set(gcf, 'PaperUnits', 'centimeters');
     set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
     saveas(gcf,[ResultFolderCell 'Cell' num2str(iCell)],'png');
    close all


end

for iCell=1:size(cellInfo,1)
    

    figure;
for iRep=1:length(RepG)
% for iRep=2:2
    I1=find(SessInfoNeed.Repetitions==RepG(iRep)&SessInfoNeed.PMTLevel==780);
    for iData=1:length(NData)
        % TempData=zscore(double(NData{iData}(iscell,:)),0,2);
        TempData=double(NData{iData}(iscell(iCell),:));
        TempData=AmpNormalize(TempData,[0 100]);

            % for iLaser=1:length(LaserG)
            for iLaser=1:8

                I2=find(SessInfoNeed.LaserPower==LaserG(iLaser));
                I3=intersect(I1,I2);
                if ~isempty(I3)
                   tempPSTH=[];

                   for iSess = 2:length(I3)
                       TempI=SessInfoNeed.PreFrame(I3(iSess)):SessInfoNeed.PostFrame(I3(iSess));
                       tempPSTH(:,iSess)=TempData(TempI);
                   end
                   PSTHLaser(:,iLaser)=squeeze(mean(tempPSTH,2));
                   % error_area(1:size(tempPSTH,1),mean(tempPSTH,2),std(tempPSTH,0,2),colorLaser(iLaser,:),0.5)
                   subplotLU(length(RepG),length(NData),iRep,iData,P);hold on


                   error_area(1:size(tempPSTH,1),mean(tempPSTH,2),ste(tempPSTH')',colorLaser(iLaser,:),0.5);
                   % if iRep==1&&iLaser==1
                   %    text(15,max(mean(tempPSTH,2)),Nlabel{iData})
                   %    % set(gca,'ylim',[0 0.5])
                   % end

                   % plot(PSTHLaser(:,iLaser),'color',colorLaser(iLaser,:));
                   % set(gca,'xtick',[0 PreImgN+0.5 PreImgN+0.5+PostImgN],'xticklabel',{'Start' num2str(PreImgN) num2str(PreImgN+PostImgN)})
                   % set(gca,'tickdir','out','clim',ClimScale(iData,:))
                   if iRep==3
                   set(gca,'xlim',[0 PreImgN+0.5+PostImgN] ,'xtick',[0 PreImgN+0.5 PreImgN+0.5+PostImgN],'xticklabel',{['-' num2str(PreImgN)] '0' num2str(PostImgN)})
                   xlabel(Nlabel{iData})
                   else
                   set(gca,'xlim',[0 PreImgN+0.5+PostImgN] ,'xtick',[0 PreImgN+0.5 PreImgN+0.5+PostImgN],'xticklabel',[])
                       
                   end

                   if iData==1
                      ylabel(['Frame #' num2str(RepG(iRep))]);
                   end
                   set(gca,'tickdir','out')
                   % colormap(ColorPN3)

                end

            end
     colormap(colorLaser)    
     b = colorbar;
     set(b,'position',[0.92 0.2 0.01 0.5]);
     ylabel(b, 'PV power');
     set(b,'ytick',[1:length(LaserG)]/length(LaserG),'yticklabel',LaserG);

    end
end
     papersizePX=[0 0 length(NData)*6 length(RepG)*5 ];
     set(gcf, 'PaperUnits', 'centimeters');
     set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
     saveas(gcf,[ResultFolderCell 'Cell' num2str(iCell) 'PMT780'],'png');
    close all


end

for iCell=1:size(cellInfo,1)
    figure;
for iRep=1:length(RepG)
% for iRep=2:2

    I1=find(SessInfoNeed.Repetitions==RepG(iRep)&SessInfoNeed.PMTLevel==780);
    for iData=1:length(NData)
        % TempData=zscore(double(NData{iData}(iscell,:)),0,2);
        TempData=double(NData{iData}(iscell(iCell),:));
        TempData=AmpNormalize(TempData,[0 100]);

            % for iLaser=1:length(LaserG)
            for iLaser=1:8

                I2=find(SessInfoNeed.LaserPower==LaserG(iLaser));
                I3=intersect(I1,I2);
                if ~isempty(I3)
                   tempPSTH=[];

                   for iSess = 2:length(I3)
                       TempI=SessInfoNeed.PreFrame(I3(iSess)):SessInfoNeed.PostFrame(I3(iSess));
                       tempPSTH(:,iSess)=TempData(TempI);
                   end
                   PSTHLaser(:,iLaser)=squeeze(mean(tempPSTH,2));
                   % error_area(1:size(tempPSTH,1),mean(tempPSTH,2),std(tempPSTH,0,2),colorLaser(iLaser,:),0.5)
                   subplotLU(length(RepG),length(NData),iRep,iData);hold on

                   error_area(1:size(tempPSTH,1),mean(tempPSTH,2),ste(tempPSTH')',colorLaser(iLaser,:),0.5);

                   % plot(PSTHLaser(:,iLaser),'color',colorLaser(iLaser,:));
                   set(gca,'xtick',[0 PreImgN+0.5 PreImgN+0.5+PostImgN],'xticklabel',{'Start' num2str(PreImgN) num2str(PreImgN+PostImgN)})
                   % set(gca,'tickdir','out','clim',ClimScale(iData,:))


                   % colormap(ColorPN3)

                end

             end
    end
end


end







