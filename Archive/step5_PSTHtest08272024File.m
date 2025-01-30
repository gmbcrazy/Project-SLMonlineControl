clear all

load('C:\Users\zhangl33\Projects\GenMatCode\Plotfun\Color\colorMapPN3.mat');

% load('C:\Users\User\Project-SLMonlineControl\subfun\Color\colorMapPN3.mat');

ProcessFolder='F:\LuSLMOnlineTest\SL0242-Ai203\08272024\SingleP\20PixelFromEdgeExc\';
DataFolder=[ProcessFolder 'Data\'];
load([ProcessFolder 'SLMIncludedIndFromIscell.mat'])
% SingP=[DataFolder 'SingleP\GPL.gpl']
SinglePSTHFolder=[DataFolder 'SinglePSTH\']
ResultFolder=SinglePSTHFolder;

mkdir(SinglePSTHFolder);

% SingPZ=[0 0 50 50 50 100 100]

BinFile=dir([DataFolder '*TSeries*GPoint*.bin'])
% % SingPZ=[0 0 0 0 50 100 100 100 100 100]
% % SingPZ=[0]
BinTable = PointLaser_files(BinFile);
TiffTable = ExpInfoMultiTiffFolder(DataFolder);

DataList = outerjoin(TiffTable, BinTable, 'Keys', 'FileID', 'MergeKeys', true);

[Pos3D,Pos3DRaw,CaData,CaDataPlane,stat,yaml]=Suite2pMultiPlaneROIToXYZ(DataFolder,2);
[cellIDMap,CellPixCount,MedCenter]=Suite2pCellIDMapFromStat(CaData.statCell,[yaml.FOVsize_PX yaml.FOVsize_PY]);
cellBoundary = CellIDMap2Boundary(cellIDMap);
yaml.Zdepth_ETL=unique(round(yaml.Zdepth_ETL));
%% In this recording file, I forgot take single image to identify the GPL Position. 
% gplFile='F:\LuSLMOnlineTest\SL0242-Ai203\08202024\Allgpl.gpl';
resultPaths = findAllFiles('F:\LuSLMOnlineTest\SL0242-Ai203\08272024\SingleP\EdgeExc\', 'GPL.gpl');
gplFile=resultPaths{1};
tbl=gpl2Table(gplFile);
% XYPos=gplXYtoPixel(tbl,yaml);
[MPPos,Z]=gplXYtoPixel(tbl,yaml);
MPName=tbl.Name;
MP3D=[MPPos yaml.Zdepth_ETL(Z(:,2))'];

%% Load Suite2p
suite2pPath = findAllFolders(DataFolder, 'suite2p');
if length(suite2pPath)==1
   CombinedSuite2p= findAllFiles([suite2pPath{1} 'combined\'], 'Fall.mat');
end
if length(CombinedSuite2p)==1
Fall=load(CombinedSuite2p{1});
cellInfo=Suite2pCellInfo(Fall);

end

PreImgN=10;
PostImgN=100;

SessInfoNeed=DataList;
LastFrame=cumsum([SessInfoNeed.totalRepetitions]);
FirstFrame=[1;LastFrame(1:end-1)+1];
clear IndStart IndEnd
for i=1:length(FirstFrame)SessInfoNeed.Laser
    IndStart(i,1)=FirstFrame(i)-PreImgN;
    IndEnd(i,1)=FirstFrame(i)+PostImgN-1;
end

% SinglePxyz=[];
% SinglePxyz=Pos3DNeed;
Laser=unique(DataList.Laser);
Point=unique(DataList.Point);
Laser(isnan(Laser))=[];
Point(isnan(Point))=[];
SinglePxyz=Pos3DNeed(Point,:);
DataList.PostSLMframe(1);


SessInfoNeed.PreFrame=IndStart;
SessInfoNeed.PostFrame=IndEnd;
SessInfoNeed(SessInfoNeed.PreFrame<0,:)=[];
NonSLMInd=find(SessInfoNeed.Laser==0&SessInfoNeed.Point>0);
SLMInd=find(SessInfoNeed.Laser>1&SessInfoNeed.Point>0);     %% Point  = 0 refers no SLM, < 0 refers to Group Stimuli. 

MarkPoint=unique(SessInfoNeed.Point(SLMInd));


colorCell = jet(length(cellBoundary));
colorCell = colorCell(randperm(length(cellBoundary)),:);


MeanImgClim=[0 1];


figure;
% Img=AmpNormalize(permute(double(CaData.PlaneMeanImg),[2 1 3]),[1 99]);
Img=AmpNormalizeDim(permute(double(CaData.PlaneMeanImg),[2 1 3]),3,[1 99]);
MultiMatrix3DPlotZ(Img,yaml.Zdepth_ETL,0.9);
caxis(MeanImgClim);
Radius=4;
colormap(gray);
set(gca,'xlim',[0 512],'ylim',[0 512],'zlim',yaml.Zdepth_ETL([1 end]),'View',[64 24],'zDir','reverse');
 % plotCellCenter3D(SinglePxyz(:,:,iPoint), Radius, [0 1 0],1.5);
% plotCellCenter3D(Pos3D, Radius, [0 1 0],1);
plotCellCenter3D(MP3D, Radius, [0 1 0],1);
labelCellCenter(MP3D, MPName,[0 1 0])
plotCellBoundary3D(cellBoundary, Pos3D(:,3),colorCell,0.5)
% labelCellCenter(Pos3D, 1:size(Pos3D,1),colorCell)
papersizePX=[0 0 16 length(yaml.Zdepth_ETL)*12];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
saveas(gcf,[ResultFolder '3DCell'],'png');
saveas(gcf,[ResultFolder '3DCell'],'fig');


CellID=1:size(Pos3D,1);

figure;
% labelCellCenter(Pos3D(I,[2 1]), CellID(I),colorCell(I,:));
MultiPlanes2DShow(Img, [], Pos3D, CellID,yaml.Zdepth_ETL, colorCell, MeanImgClim)
papersizePX=[0 0 length(yaml.Zdepth_ETL)*11 10];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
saveas(gcf,[ResultFolder '2DCell'],'png');
saveas(gcf,[ResultFolder '2DCell'],'fig');




SLMCellID=[18 22 19 25 3 44 47 -50 65 -66 -64]  %Level4

LaserG=Laser
% RepG=unique(SessInfoNeed.Repetitions)
% PmtG=unique(SessInfoNeed.PMTLevel)

deltaFoF=F2deltaFoF(Fall.F,Fall.Fneu,Fall.ops.fs);


NData={deltaFoF', Fall.spks};
Nlabel={'DeltaF', 'Spks'}
planeG=unique(cellInfo.iplane);
for iplane=1:length(planeG)
    PlaneC(iplane)=max(find(cellInfo.iplane==planeG(iplane)));
end

iscell=find(Fall.iscell(:,1)==1);
ClimScale=[-5 5;-40 40]
ClimScale=[-1 1;-1 1]
ClimScale=[-2 2;-5 5]
ClimScale=[-0.3 0.3;-0.3 0.3]

ClimScale=[-5 5;-15 15]
   P.xLeft=0.1;        %%%%%%Left Margin
   P.xRight=0.1;       %%%%%%Right Margin
   P.yTop=0.02;         %%%%%%Top Margin
   P.yBottom=0.1;      %%%%%%Bottom Margin
   P.xInt=0.01;         %%%%%%Width-interval between subplots
   P.yInt=0.03;         %%%%%%Height-interval between subplots

SLMframe=[50];


    % I1=NonSLMInd;
    figure;
    for iData=1:length(NData)
        % TempData=zscore(double(NData{iData}(iscell,:)),0,2);
        TempData=double(NData{iData}(iscell,:));
        % for iCell=1:size(TempData)
        %     TempData(iCell,:)=AmpNormalize(TempData(iCell,:),[0 100]);
        % end
        % 

        % for iCell=1:size(cellInfo,1)
            for iLaser=1:length(LaserG)
               
                I2=find(SessInfoNeed.Laser==LaserG(iLaser));
                I3=I2;
                if ~isempty(I3)
                   tempPSTH=[];

                   for iSess = 1:length(I3)
                       TempI=SessInfoNeed.PreFrame(I3(iSess)):SessInfoNeed.PostFrame(I3(iSess));
                       tempPSTH(:,:,iSess)=TempData(:,TempI);
                   end
                   tempPSTH=squeeze(mean(tempPSTH,3));
                   BaseLine=repmat(mean(tempPSTH(:,1:PreImgN),2),1,size(tempPSTH,2));
                   tempPSTH=(tempPSTH-BaseLine)./BaseLine;



                   subplotLU(length(NData),length(LaserG),iData,iLaser,P);
                   imagesc(tempPSTH);hold on;
                   if iData==2
                   set(gca,'xlim',[0 PreImgN+0.5+PostImgN] ,'xtick',SLMframe+PreImgN+0.5)
                   xlabel(['PV-Power ' num2str(LaserG(iLaser))])
                   else
                   set(gca,'xlim',[0 PreImgN+0.5+PostImgN] ,'xtick',SLMframe+PreImgN+0.5)
                       
                   end
                   set(gca,'tickdir','out','clim',ClimScale(iData,:))
                   colormap(ColorPN3)
                   if iLaser >1 
                      set(gca,'ytick',[]);
                   elseif   iLaser == 1 
                      ylabel('CellsID');
                      set(gca,'ytick',PlaneC+0.5,'yticklabel',{'Plane1' 'Plane2' 'Plane3'})

                   else
                   end

                   for iplane=1:length(PlaneC)
                       plot([0 PreImgN+0.5+PostImgN],zeros(1,2)+PlaneC(iplane)+0.5,'g-')
                   end
                   % plot([PreImgN PreImgN]+0.5,[0 length(iscell)],'k-')

                end

            end
         % b = colorbar;
         % set(b,'position',[0.97 0.7-(iData-1)*0.5 0.002 0.2]);
         % ylabel(b, Nlabel{iData});

        % end
    end
    papersizePX=[0 0 length(LaserG)*5 length(NData)*8];
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
    saveas(gcf,[ResultFolder 'AllNeuroSessReptition1Trial'],'png');
    saveas(gcf,[ResultFolder 'AllNeuroSessReptition1Trial'],'fig');

    close all






PlaneNumStart=[1 PlaneC(1:end-1)+1]
jet=colormap("jet");
colorLaser=jet(1:size(jet,1)/length(LaserG):size(jet,1),:);
close all

% I1=intersect(NonSLMInd,SingleRep);
figure;
for iData=1:length(NData)
        % TempData=zscore(double(NData{iData}(iscell,:)),0,2);
        for iPlane = 1:length(PlaneC)


        TempData=double(NData{iData}(iscell(PlaneNumStart(iPlane):PlaneC(iPlane)),:));
        for iCell=1:size(TempData)
            TempData(iCell,:)=AmpNormalize(TempData(iCell,:),[0 100]);
        end

        subplotLU(length(NData),length(PlaneC),iData,iPlane,P);

        % for iCell=1:size(cellInfo,1)
            for iLaser=1:length(LaserG)
               
                I2=find(SessInfoNeed.Laser==LaserG(iLaser));
                I3=I2;
                if ~isempty(I3)
                   tempPSTH=[];
                   for iSess = 2:length(I3)
                       TempI=SessInfoNeed.PreFrame(I3(iSess)):SessInfoNeed.PostFrame(I3(iSess));
                       tempPSTH(:,:,iSess)=TempData(:,TempI);
                   end

                   % for iCell=1:size(tempPSTH,1)
                   %     tempPSTH(iCell,:,:)=AmpNormalize(tempPSTH(iCell,:,:),[0 100]);
                   % end
                   tempPSTH=squeeze(mean(tempPSTH,3));
                   BaseLine=repmat(mean(tempPSTH(:,1:PreImgN),2),1,size(tempPSTH,2));
                   tempPSTH=tempPSTH-BaseLine;

                   error_area(1:size(tempPSTH,2),mean(tempPSTH,1),ste(tempPSTH),colorLaser(iLaser,:),0.1);

                   % imagesc(tempPSTH);hold on;
                   if iData==2
                   set(gca,'xlim',[0 PreImgN+0.5+PostImgN] ,'xtick',[0 PreImgN+0.5 PreImgN+0.5+PostImgN],'xticklabel',{['-' num2str(PreImgN)] '0' num2str(PostImgN)})
                   xlabel(['Plane ' num2str(iPlane)])
                   else
                   set(gca,'xlim',[0 PreImgN+0.5+PostImgN] ,'xtick',[0 PreImgN+0.5 PreImgN+0.5+PostImgN],'xticklabel',[])

                   end
                   set(gca,'tickdir','out')
                   set(gca,'ylim',[-0.01 0.02],'ytick',[-0.02:0.02:0.04])
                   % colormap(ColorPN3)
                   % if iLaser >1 
                   %    set(gca,'ytick',[]);
                   % elseif   iLaser == 1 
                   %    ylabel('CellsID');
                   %    set(gca,'ytick',PlaneC+0.5,'yticklabel',{'Plane1' 'Plane2' 'Plane3'})
                   % 
                   % else
                   % end
                   % 
                   % for iplane=1:length(PlaneC)
                   %     plot([0 PreImgN+0.5+PostImgN],zeros(1,2)+PlaneC(iplane)+0.5,'g-')
                   % end
                   % plot([PreImgN PreImgN]+0.5,[0 0.05],'k-')

                end

            end
        
        
        if iPlane==1
           ylabel(Nlabel{iData})
        end
        
        
        end
         % colormap(colorLaser)
     % set(b,'tick',[1:length(LaserG)]/length(LaserG),'ticklabel',LaserG);

        % end
end

     colormap(colorLaser)    
     b = colorbar;
     set(b,'position',[0.92 0.2 0.01 0.5]);
     ylabel(b, 'PV power');
     set(b,'ytick',[1:length(LaserG)]/length(LaserG),'yticklabel',LaserG);
     papersizePX=[0 0 length(PlaneC)*6 length(NData)*5];
     set(gcf, 'PaperUnits', 'centimeters');
     set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
     saveas(gcf,[ResultFolder 'AveAllNeuroSessReptition1Trial'],'png');
    close all



for iData=1:length(NData)
figure;
hold on;
    for iCell=1:size(cellInfo,1)

        % TempData=zscore(double(NData{iData}(iscell,:)),0,2);
        TempData=double(NData{iData}(iscell(iCell),:));
        TempData=AmpNormalize(TempData,[0 100]);
        % plot(TempData-iCell*1,'color',[0.2 0.2 0.2]);
        plot(smooth(TempData,3)-iCell*1,'color',[0.2 0.2 0.2]);

        EndPlane=ismember(iCell,PlaneC)
        if EndPlane
           plot([0 length(TempData)],zeros(1,2)-iCell,'g-');
        end

    end
    for iS=1:size(SessInfoNeed,1)
       [~,Laserj]=ismember(SessInfoNeed.Laser(iS),LaserG);
       % plot(zeros(1,2)+SessInfoNeed.PreFrame(iS)+PreImgN,[-iCell 0],'color',colorLaser(Laserj,:));



         SessStart=SessInfoNeed.PreFrame(iS)+PreImgN;
         RepetitionFameN=80;
         for iRep = 1:length(SLMframe)
             S1=SessStart+SLMframe(iRep);
             if iRep==1
                plot(zeros(1,2)+S1,[-iCell 0],':','color',colorLaser(Laserj,:));
             else
                plot(zeros(1,2)+S1,[-iCell 0],'color',colorLaser(Laserj,:));
             end
         end
       % text(SessInfoNeed.PreFrame(iS)+PreImgN,0,['Laser' num2str(LaserG(Laserj))],'HorizontalAlignment','center');
    end

    set(gca,'ytick',[-length(iscell):1:-1]+0.5,'yticklabel',abs([-length(iscell):1:-1]))
    ylabel('Cell IDs')
     colormap(colorLaser)    
     b = colorbar;
     set(b,'position',[0.9 0.2 0.03 0.5]);
     ylabel(b, 'PV power');
     set(b,'ytick',[1:length(LaserG)]/length(LaserG),'yticklabel',LaserG);
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



%% Ave SLM targets illustration
close all
ResultFolderCell=[ResultFolder 'SLMtarget\'];
mkdir(ResultFolderCell);
% for iCell=1:size(cellInfo,1)

   P.xLeft=0.1;        %%%%%%Left Margin
   P.xRight=0.1;       %%%%%%Right Margin
   P.yTop=0.02;         %%%%%%Top Margin
   P.yBottom=0.1;      %%%%%%Bottom Margin
   P.xInt=0.01;         %%%%%%Width-interval between subplots
   P.yInt=0.03;         %%%%%%Height-interval between subplots
% I1=intersect(NonSLMInd,SingleRep);
close all
for jCell=1:length(SLMCellID)

    iCell=SLMCellID(jCell);
    if iCell<0
       continue;
    end
    figure;
% for iRep=2:2
    for iData=1:length(NData)
        % TempData=zscore(double(NData{iData}(iscell,:)),0,2);
        TempData=double(NData{iData}(iscell(iCell),:));
        TempData=AmpNormalize(TempData,[0 100]);

            for iLaser=1:length(LaserG)
            % for iLaser=1:8

                if LaserG(iLaser)==0
                   I2=find(SessInfoNeed.Laser==LaserG(iLaser));      
                else
                   I2=find(SessInfoNeed.Laser==LaserG(iLaser)&SessInfoNeed.Point==MarkPoint(jCell));

                end

                I3=I2;
                if ~isempty(I3)
                   tempPSTH=[];

                   for iSess = 1:length(I3)
                       TempI=SessInfoNeed.PreFrame(I3(iSess)):SessInfoNeed.PostFrame(I3(iSess));
                       tempPSTH(:,iSess)=TempData(TempI);
                   end
                   PSTHLaser(:,iLaser)=squeeze(mean(tempPSTH,2));
                   % error_area(1:size(tempPSTH,1),mean(tempPSTH,2),std(tempPSTH,0,2),colorLaser(iLaser,:),0.5)
                   subplotLU(1,length(NData),1,iData,P);hold on
                   BaseLine=repmat(mean(tempPSTH(1:PreImgN,:),1),size(tempPSTH,1),1);
                   tempPSTH=tempPSTH-BaseLine;


                   error_area(1:size(tempPSTH,1),mean(tempPSTH,2),ste(tempPSTH')',colorLaser(iLaser,:),0.4,'-',0.5);
                   % if iRep==1&&iLaser==1
                   %    text(15,max(mean(tempPSTH,2)),Nlabel{iData})
                   %    % set(gca,'ylim',[0 0.5])
                   % end
                   set(gca,'ylim',[-0.02 0.12]);
                   hold on;
                   for jf=1:length(SLMframe)
                   plot(PreImgN+[SLMframe(jf) SLMframe(jf)],[-0.1 0.3],'k:');
                   end


                   % plot(PSTHLaser(:,iLaser),'color',colorLaser(iLaser,:));
                   % set(gca,'xtick',[0 PreImgN+0.5 PreImgN+0.5+PostImgN],'xticklabel',{'Start' num2str(PreImgN) num2str(PreImgN+PostImgN)})
                   % set(gca,'tickdir','out','clim',ClimScale(iData,:))
                   % if iRep==3
                   set(gca,'xlim',[0 PreImgN+0.5+PostImgN] ,'xtick',[0 PreImgN+0.5 PreImgN+0.5+PostImgN],'xticklabel',{['-' num2str(PreImgN)] '0' num2str(PostImgN)})
                   xlabel(Nlabel{iData})
                   % else
                   % set(gca,'xlim',[0 PreImgN+0.5+PostImgN] ,'xtick',[0 PreImgN+0.5 PreImgN+0.5+PostImgN],'xticklabel',[])
                   % 
                   % end
                   set(gca,'tickdir','out')
                   set(gca,'ylim',[-0.2 0.6],'ytick',[-0.2:0.2:0.6])
                   % if iData==1
                   %    ylabel(['Frame #' num2str(RepG(iRep))]);
                   % end

                   % colormap(ColorPN3)

                end

            end
     colormap(colorLaser)    
     b = colorbar;
     set(b,'position',[0.92 0.2 0.01 0.5]);
     ylabel(b, 'PV power');
     set(b,'ytick',[1:length(LaserG)]/length(LaserG),'yticklabel',LaserG);

    end
     papersizePX=[0 0 length(NData)*8 5 ];
     set(gcf, 'PaperUnits', 'centimeters');
     set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
     saveas(gcf,[ResultFolderCell 'SLMTarget' num2str(jCell)],'png');
     % close all


end


%% All Cells
jCell=1;
ResultFolderCell=[ResultFolder 'Cells\'];
mkdir(ResultFolderCell);
close all
for iCell=1:size(cellInfo,1)

        figure;
 
% for iRep=2:2
    for iData=1:length(NData)
        % TempData=zscore(double(NData{iData}(iscell,:)),0,2);
        TempData=double(NData{iData}(iscell(iCell),:));
        TempData=AmpNormalize(TempData,[0 100]);
        for jCell=1:length(SLMCellID)

            for iLaser=1:length(LaserG)
            % for iLaser=1:8

                if LaserG(iLaser)==0
                   I2=find(SessInfoNeed.Laser==LaserG(iLaser));      
                else
                   I2=find(SessInfoNeed.Laser==LaserG(iLaser)&SessInfoNeed.Point==MarkPoint(jCell));
                end

                I3=I2;
                if ~isempty(I3)
                   tempPSTH=[];

                   for iSess = 1:length(I3)
                       TempI=SessInfoNeed.PreFrame(I3(iSess)):SessInfoNeed.PostFrame(I3(iSess));
                       tempPSTH(:,iSess)=TempData(TempI);
                   end
                   PSTHLaser(:,iLaser)=squeeze(mean(tempPSTH,2));
                   % error_area(1:size(tempPSTH,1),mean(tempPSTH,2),std(tempPSTH,0,2),colorLaser(iLaser,:),0.5)
                   subplotLU(length(SLMCellID),length(NData),jCell,iData,P);hold on
                   BaseLine=repmat(mean(tempPSTH(1:PreImgN,:),1),size(tempPSTH,1),1);
                   tempPSTH=tempPSTH-BaseLine;


                   error_area(1:size(tempPSTH,1),mean(tempPSTH,2),ste(tempPSTH')',colorLaser(iLaser,:),0.4,'-',0.5);
                   % if iRep==1&&iLaser==1
                   %    text(15,max(mean(tempPSTH,2)),Nlabel{iData})
                   %    % set(gca,'ylim',[0 0.5])
                   % end
                   set(gca,'ylim',[-0.02 0.12]);
                   hold on;
                   for jf=1:length(SLMframe)
                   plot(PreImgN+[SLMframe(jf) SLMframe(jf)],[-0.1 0.3],'k:');
                   end


                   % plot(PSTHLaser(:,iLaser),'color',colorLaser(iLaser,:));
                   % set(gca,'xtick',[0 PreImgN+0.5 PreImgN+0.5+PostImgN],'xticklabel',{'Start' num2str(PreImgN) num2str(PreImgN+PostImgN)})
                   % set(gca,'tickdir','out','clim',ClimScale(iData,:))
                   if jCell==length(SLMCellID)
                   set(gca,'xlim',[0 PreImgN+0.5+PostImgN] ,'xtick',[0 PreImgN+0.5 PreImgN+0.5+PostImgN],'xticklabel',{['-' num2str(PreImgN)] '0' num2str(PostImgN)})
                   xlabel(Nlabel{iData})
                   else
                   set(gca,'xlim',[0 PreImgN+0.5+PostImgN] ,'xtick',[0 PreImgN+0.5 PreImgN+0.5+PostImgN],'xticklabel',[])

                   end
                   set(gca,'tickdir','out')
                   set(gca,'ylim',[-0.1 0.3],'ytick',[-0.1:0.1:0.3])
                   % if iData==1
                   %    ylabel(['Frame #' num2str(RepG(iRep))]);
                   % end

                   % colormap(ColorPN3)

                end

            end
     
        end
            colormap(colorLaser)    
     b = colorbar;
     set(b,'position',[0.92 0.2 0.01 0.5]);
     ylabel(b, 'PV power');
     set(b,'ytick',[1:length(LaserG)]/length(LaserG),'yticklabel',LaserG);

    end
     papersizePX=[0 0 length(NData)*8 length(SLMCellID)*5];
     set(gcf, 'PaperUnits', 'centimeters');
     set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
     saveas(gcf,[ResultFolderCell 'Cell' num2str(iCell)],'png');
     close all


end


















nPlane=3
PreInd=40:60
PostInd=62:64;
PlaneZ=confSet.ETL+confSet.scan_Z
% [0 50 100];
MeanImgClim=[-150 150];








for iLaser=1:length(Laser)
    Ind1=find(LaserFile==Laser(iLaser));
             % figure;
    tic
    for iPoint = 1:length(Point)
         Ind2=find(PointFile==Point(iPoint));
         IndNeed=intersect(Ind1,Ind2);
         BinFolderList=BinFileFolder(IndNeed);

         for iTrial=1:length(BinFolderList)
             Temp1(:,:,:,iTrial)=MeanFrameIndMultiTiffs([DataFolder BinFolderList{iTrial} '\'],nPlane,PreInd);
             Temp2(:,:,:,iTrial)=MeanFrameIndMultiTiffs([DataFolder BinFolderList{iTrial} '\'],nPlane,PostInd);

             if iTrial==1
             xmlFile=dir([DataFolder BinFolderList{iTrial} '\*MarkPoints.xml']);
             xmlInfo=xml2struct([xmlFile(1).folder '\' xmlFile(1).name]);
             % pointInfo=xmlInfo.PVMarkPointSeriesElements.PVMarkPointElement.PVGalvoPointElement.Point.Attributes;
             pointInfo=xmlInfo.PVMarkPointSeriesElements.PVMarkPointElement.PVGalvoPointElement.Point{1}.Attributes;

             X=str2num(pointInfo.X)*512;
             Y=str2num(pointInfo.Y)*512;
             SinglePxyz(iPoint,:)=[Y X SingPZ(iPoint)];

             end
         end
         Temp3=mean(Temp2,4)-mean(Temp1,4);
         Temp3=SmoothDecDim3(Temp3,1.5);
    
         % subplot(3,3,iPoint)
%           figure;
%          MultiMatrix3DPlotZ(Temp3,PlaneZ,0.9);
%         caxis(MeanImgClim);
%          Radius=8;
%         colormap(ColorPN3);
% % set(gca,'xlim',[0 512],'ylim',[0 512],'zlim',[-80 80],'clim',ClimImg*10,'yDir','reverse');
%         set(gca,'xlim',[0 512],'ylim',[0 512],'zlim',PlaneZ([1 end]),'View',[64 24],'zDir','reverse');
%           plotCellCenter3D(SinglePxyz(iPoint,:), Radius, [0 1 0],2.5);

         PSTH{iLaser,iPoint}=Temp3;
    end
    toc
end

close all
papersizePX=[0 0 22 18];

for iPoint = 1:length(Point)
    figure;

     for iLaser=1:length(Laser)
         subplot(1,length(Laser),iLaser)
         MultiMatrix3DPlotZ(PSTH{iLaser,iPoint},PlaneZ,0.9);
         caxis(MeanImgClim);
         Radius=10;
       colormap(ColorPN3);
         set(gca,'xlim',[0 512],'ylim',[0 512],'zlim',PlaneZ([1 end]),'View',[64 24],'zDir','reverse');
         plotCellCenter3D(SinglePxyz(:,:,iPoint), Radius, [0 1 0],1.5);

     
     end
%       set(gcf, 'PaperUnits', 'centimeters');
%       set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
%       saveas(gcf,[SinglePSTHFolder 'Point' num2str(Point(iPoint))],'png'); 
%       saveas(gcf,[SinglePSTHFolder 'Point' num2str(Point(iPoint))],'fig'); 

     % close all
end


PreInd=45:60
PostInd=62:63;
MeanImgClim=[-300 300]
TiffFolderName='TSeries-04222024-0926-';
IndexFolderName=[102:105 107:109];
clear Temp1 Temp2
for iTrial=1:length(IndexFolderName)
    tic
    tempFile= [DataFolder TiffFolderName num2str(IndexFolderName(iTrial)) '\'];
    Temp1(:,:,:,iTrial)=MeanFrameIndMultiTiffs(tempFile,nPlane,PreInd);
    Temp2(:,:,:,iTrial)=MeanFrameIndMultiTiffs(tempFile,nPlane,PostInd);
    toc
end
 Temp3=mean(Temp2,4)-mean(Temp1,4);
 Temp3=SmoothDecDim3(Temp3,1.5);

figure;
subplot('position',[0.01 0.01 0.88 0.95])
MultiMatrix3DPlotZ(Temp3,PlaneZ,0.9);
caxis(MeanImgClim);
Radius=8;
colormap(ColorPN3);
% SinglePxyz=squeeze(SinglePxyz);
% SinglePxyz=SinglePxyz';
% set(gca,'xlim',[0 512],'ylim',[0 512],'zlim',[-80 80],'clim',ClimImg*10,'yDir','reverse');
set(gca,'xlim',[0 512],'ylim',[0 512],'zlim',PlaneZ([1 end]),'View',[64 24],'zDir','reverse');
plotCellCenter3D(SinglePxyz, Radius, [0 1 0],2.5);
% 




% % nPlane=3;
% % for iPlane=1:nPlane
% % 
% % for iframe=1:length(TiffFile)
% %     X(:,:,iframe,iPlane)=imread([TiffFile(iframe).folder '\' TiffFile(iframe).name],iPlane);
% % end
% %     MeanImage(:,:,iPlane) = mean(X(:,:,:,iPlane),3);
% % end



PreInd=15:25
PostInd=26:28;

% PreInd=2:11;
% PauseInd=11;
% PostInd=13:26;
% PlaneZ=[0 50 100];
clear PSTHtemp;
% BinFile=dir([DataFolder '*TSeries*Laser1.12Point6.bin'])
FrameTotal=(80)*nPlane
% Skip=PauseInd*3+1;




Ly=512;
Lx=512;
% PreInd=40:60
% PostInd=61:63;


for iLaser=1:length(Laser)
    Ind1=find(LaserFile==Laser(iLaser));
             % figure;
    tic
    for iPoint = 1:length(Point)
         Ind2=find(PointFile==Point(iPoint));
         IndNeed=intersect(Ind1,Ind2);
         BinFolderList=BinFileFolder(IndNeed);
         clear PSTHtemp;
         for iTrial=1:length(BinFolderList)
             iFile=IndNeed(iTrial);
             fileID=fopen([BinFile(iFile).folder '\' BinFile(iFile).name]);
             clear Y

            BinData=fread(fileID,'uint16');

            BinDataAll=BinData(1:Ly*Lx*FrameTotal);   
            % FileCount(2)=FileCount(2)+1;


             % FrameTotal=floor(length(BinDataAll)/Ly/Lx);
            fclose(fileID);

            XBin=reshape(BinDataAll,Ly,Lx,FrameTotal);

% clear Y
            Y=zeros(Lx,Ly,floor(FrameTotal/nPlane),nPlane);
            for i=1:nPlane
                Y(:,:,1:FrameTotal/nPlane,i)=XBin(1:Lx,1:Ly,i:nPlane:size(XBin,3));
            end

            Temp11(:,:,:,iTrial)=mean(Y(:,:,PreInd,:),3);
            Temp22(:,:,:,iTrial)=mean(Y(:,:,PostInd-1,:),3);
            PSTHtemp(:,:,:,iTrial)=mean(Y(:,:,PostInd-1,:),3)-mean(Y(:,:,PreInd,:),3);

         end
         PSTHtemp=squeeze(mean(PSTHtemp,4));
         PSTHtemp=permute(PSTHtemp,[2,1,3]);
         PSTHtemp=SmoothDecDim3(PSTHtemp,1.5);

         PSTHBin{iLaser,iPoint}=PSTHtemp;
    end
    toc
end


for iLaser=1:length(Laser)
    Ind1=find(LaserFile==Laser(iLaser));
             % figure;
    tic
    for iPoint = 1:length(Point)
         Ind2=find(PointFile==Point(iPoint));
         IndNeed=intersect(Ind1,Ind2);
         BinFolderList=BinFileFolder(IndNeed);
         clear PSTHtemp;
         for iTrial=1:length(BinFolderList)
             iFile=IndNeed(iTrial);
             fileID=fopen([BinFile(iFile).folder '\' BinFile(iFile).name]);
             clear Y

             PreData = Suite2pSingleChBin2Frame([BinFile(iFile).folder '\' BinFile(iFile).name], Ly, Lx, nPlane, PreInd);
             PostData= Suite2pSingleChBin2Frame([BinFile(iFile).folder '\' BinFile(iFile).name], Ly, Lx, nPlane, PostInd-1);

            % Temp11(:,:,:,iTrial)=mean(PreData,3);
            % Temp22(:,:,:,iTrial)=mean(PostData,3);
            PSTHtemp(:,:,:,iTrial)=mean(PostData,3)-mean(PreData,3);
         end

         PSTHtemp=squeeze(mean(PSTHtemp,4));
         % PSTHtemp=permute(PSTHtemp,[2,1,3]);
         PSTHtemp=SmoothDecDim3(PSTHtemp,1.5);

         PSTHBin{iLaser,iPoint}=PSTHtemp;
    end
    toc
end




for iTrial=1:5
figure;
subplot(1,2,1);imagesc(Temp2(:,:,1,iTrial))
subplot(1,2,2);imagesc(Temp22(:,:,1,iTrial)')
end


close all
papersizePX=[0 0 22 18*length(Laser)];
for iPoint = 1:length(Point)
    figure;
     for iLaser=1:length(Laser)
         % subplot(length(Laser),2,(iLaser-1)*2+1)
         % MultiMatrix3DPlotZ(PSTH{iLaser,iPoint},PlaneZ,0.9);
         % caxis(MeanImgClim);
         % Radius=10;
         % colormap(ColorPN3);
         % set(gca,'xlim',[0 512],'ylim',[0 512],'zlim',PlaneZ([1 end]),'View',[64 24],'zDir','reverse');
         % plotCellCenter3D(SinglePxyz(iPoint,:), Radius, [0 1 0],1.5);
         % ylabel('Tiff Data')

         subplot(length(Laser),2,(iLaser-1)*2+2)
         MultiMatrix3DPlotZ(PSTHBin{iLaser,iPoint},PlaneZ,0.9);
         caxis(MeanImgClim);
         Radius=10;
         colormap(ColorPN3);
         set(gca,'xlim',[0 512],'ylim',[0 512],'zlim',PlaneZ([1 end]),'View',[64 24],'zDir','reverse');
         % plotCellCenter3D(SinglePxyz(iPoint,:), Radius, [0 1 0],1.5);
         plotCellCenter3D(SinglePxyz(iPoint,:), Radius, [0 1 0],1.5);

         ylabel('Bin Data')

     end
      set(gcf, 'PaperUnits', 'centimeters');
      set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
      % saveas(gcf,[SinglePSTHFolder 'Point' num2str(Point(iPoint))],'png'); 
      % saveas(gcf,[SinglePSTHFolder 'Point' num2str(Point(iPoint))],'fig'); 

     % close all
end



close all
papersizePX=[0 0 22*2 18*length(Laser)];
for iPoint = 1:length(Point)
    figure;
     for iLaser=1:length(Laser)
         % subplot(length(Laser),2,(iLaser-1)*2+1)
         % MultiMatrix3DPlotZ(PSTH{iLaser,iPoint},PlaneZ,0.9);
         % caxis(MeanImgClim);
         % Radius=10;
         % colormap(ColorPN3);
         % set(gca,'xlim',[0 512],'ylim',[0 512],'zlim',PlaneZ([1 end]),'View',[64 24],'zDir','reverse');
         % plotCellCenter3D(SinglePxyz(:,:,iPoint), Radius, [0 1 0],1.5);
         % ylabel('Tiff Data')

         subplot(1,length(Laser),iLaser)
         MultiMatrix3DPlotZ(PSTHBin{iLaser,iPoint},PlaneZ,0.9);
         caxis(MeanImgClim);
         Radius=10;
         colormap(ColorPN3);
         set(gca,'xlim',[0 512],'ylim',[0 512],'zlim',PlaneZ([1 end]),'View',[64 24],'zDir','reverse');
         plotCellCenter3D(SinglePxyz(:,:,iPoint), Radius, [0 1 0],1.5);
         ylabel('Bin Data')

     end
      set(gcf, 'PaperUnits', 'centimeters');
      set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
      saveas(gcf,[SinglePSTHFolder 'Point' num2str(Point(iPoint))],'png'); 
      saveas(gcf,[SinglePSTHFolder 'Point' num2str(Point(iPoint))],'fig'); 

     close all
end






