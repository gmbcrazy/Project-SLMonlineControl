
%% In case of two channel recording, require to move Ch2 tiff to another folder for easier processing.
% clear all
% 
% path='D:\LuSLMOnlineTest\SL0340\05302024\';
% DestiFolder='D:\LuSLMOnlineTest\SL0340\05302024\Ch2\';
% Folder='TSeries';
% NeedFile='tif';
% keyWord='Ch2';
% copyFormat(path,Folder,DestiFolder,NeedFile,keyWord)
% 
% 
% Folder='TSeries';
% NeedFile='xml';
% keyWord='TSeries';
% copyFormat(path,Folder,DestiFolder,NeedFile,keyWord)
% NeedFile='env';
% keyWord='TSeries';
% copyFormat(path,Folder,DestiFolder,NeedFile,keyWord)
% NeedFile='ome';
% keyWord='TSeries';
% copyFormat(path,Folder,DestiFolder,NeedFile,keyWord)
% 




%% 

clear all
ProcessFolder='F:\LuSLMOnlineTest\SL0242-Ai203\08192024\';
PreImgN=10;   %%Define the frame Nums before each imaging file start to show, coming from the previous file
PostImgN=90;  %%Define the frame Nums after each imaging file start to show
ResultFolder=[ProcessFolder 'Result\'];
mkdir(ResultFolder)
SLMframe=[0 1 51];


load('C:\Users\zhangl33\Projects\GenMatCode\Plotfun\Color\colorMapPN3.mat');

[Pos3D,Pos3DRaw,CaData,CaDataPlane,stat,yaml]=Suite2pMultiPlaneROIToXYZ(ProcessFolder);
[cellIDMap,CellPixCount,MedCenter]=Suite2pCellIDMapFromStat(CaData.statCell,[yaml.FOVsize_PX yaml.FOVsize_PY]);
cellBoundary = CellIDMap2Boundary(cellIDMap);
yaml.Zdepth_ETL=unique(round(yaml.Zdepth_ETL))




%% Load Suite2p
suite2pPath = findAllFolders(ProcessFolder, 'suite2p');
if length(suite2pPath)==1
   CombinedSuite2p= findAllFiles([suite2pPath{1} 'combined\'], 'Fall.mat');
end
if length(CombinedSuite2p)==1
Fall=load(CombinedSuite2p{1});
cellInfo=Suite2pCellInfo(Fall);
end

%% Find Tiff Folder
TiffFolder=dir([ProcessFolder 'TSeries*'])
for i=1:size(TiffFolder,1)
    Session(i)=str2num(TiffFolder(i,1).name(end-2:end));
    SessionNum(i)=length(dir([TiffFolder(i,1).folder '\' TiffFolder(i,1).name '\*.tif']))/3;

end


SessInfoPath=ProcessFolder;
SessInfo=readtable([SessInfoPath 'SessionInfo.xlsx'])
[~,I1,I2]=intersect(Session,SessInfo.Session);
SessInfoNeed=SessInfo(I2,:);
LastFrame=cumsum([SessInfoNeed.TotalRepetition]);
FirstFrame=[1;LastFrame(1:end-1)+1];
% % LastFrame=7500+cumsum([SessInfoNeed.TotalRepetition]);
% % FirstFrame=[7501;LastFrame(1:end-1)+1];



clear IndStart IndEnd
for i=1:length(FirstFrame)
    IndStart(i,1)=FirstFrame(i)-PreImgN;
    IndEnd(i,1)=FirstFrame(i)+PostImgN-1;
end



SessInfoNeed.PreFrame=IndStart;
SessInfoNeed.PostFrame=IndEnd;
SessInfoNeed(SessInfoNeed.PreFrame<0,:)=[];
NonSLMInd=find(SessInfoNeed.Power==0&SessInfoNeed.Point>0);
SLMInd=find(SessInfoNeed.Power>220&SessInfoNeed.Point>0);     %% Point  = 0 refers no SLM, < 0 refers to Group Stimuli. 

MarkPoint=unique(SessInfoNeed.Point(SLMInd));

% SingleRep=find(SessInfoNeed.RepeatTimes==1);
% MulRep=find(SessInfoNeed.RepeatTimes>1);
colorCell = jet(length(cellBoundary));
colorCell = colorCell(randperm(length(cellBoundary)),:);


 MeanImgClim=[0 1]
 figure;
 Img=AmpNormalize(permute(double(CaData.PlaneMeanImg),[2 1 3]),[1 99.5]);
 MultiMatrix3DPlotZ(Img,yaml.Zdepth_ETL,0.9);
 caxis(MeanImgClim);
 Radius=4;
 colormap(gray);
 set(gca,'xlim',[0 512],'ylim',[0 512],'zlim',yaml.Zdepth_ETL([1 end]),'View',[64 24],'zDir','reverse');
 % plotCellCenter3D(SinglePxyz(:,:,iPoint), Radius, [0 1 0],1.5);
% plotCellCenter3D(Pos3D, Radius, [0 1 0],1);


plotCellBoundary3D(cellBoundary, Pos3D(:,3),colorCell)
% labelCellCenter(Pos3D, 1:size(Pos3D,1),colorCell)
    papersizePX=[0 0 16 length(yaml.Zdepth_ETL)*12];
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
    saveas(gcf,[ResultFolder '3DCell'],'png');
    saveas(gcf,[ResultFolder '3DCell'],'fig');



figure;
for iplane=1:length(yaml.Zdepth_ETL)
subplot(1,3,iplane)
imagesc(Img(:,:,iplane)');
I=find(abs(Pos3D(:,3)-yaml.Zdepth_ETL(iplane))<0.1);
CellID=1:size(Pos3D,1);
% plotCellBoundary(cellBoundary(I), CellID(I))
colormap(gray);
caxis(MeanImgClim);

labelCellCenter(Pos3D(I,[2 1]), CellID(I),colorCell(I,:));
end 
% plotCellBoundary3D(cellBoundary, Pos3D(:,3))
% labelCellCenter(Pos3D, 1:size(Pos3D,1))
    papersizePX=[0 0 length(yaml.Zdepth_ETL)*11 10];
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
    saveas(gcf,[ResultFolder '2DCell'],'png');
    saveas(gcf,[ResultFolder '2DCell'],'fig');



%%Last one:313 461 450
% SLMCellID=[15 66 70 69]  %Level1
% SLMCellID=[1 74 91 60]  %Level2
% SLMCellID=[27 225 161 166]  %Level3
%P3 42 40
SLMCellID=[27 18 42 47 51 87 81 80 79]  %Level4

LaserG=unique(SessInfoNeed.Power)
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
               
                I2=find(SessInfoNeed.Power==LaserG(iLaser));
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
               
                I2=find(SessInfoNeed.Power==LaserG(iLaser));
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
       [~,Laserj]=ismember(SessInfoNeed.Power(iS),LaserG);
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
    figure;
% for iRep=2:2
    for iData=1:length(NData)
        % TempData=zscore(double(NData{iData}(iscell,:)),0,2);
        TempData=double(NData{iData}(iscell(iCell),:));
        TempData=AmpNormalize(TempData,[0 100]);

            for iLaser=1:length(LaserG)
            % for iLaser=1:8

                if LaserG(iLaser)==0
                   I2=find(SessInfoNeed.Power==LaserG(iLaser));      
                else
                   I2=find(SessInfoNeed.Power==LaserG(iLaser)&SessInfoNeed.Point==MarkPoint(jCell));

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
                   I2=find(SessInfoNeed.Power==LaserG(iLaser));      
                else
                   I2=find(SessInfoNeed.Power==LaserG(iLaser)&SessInfoNeed.Point==MarkPoint(jCell));
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




for iCell=1:size(cellInfo,1)


    figure;
% for iRep=2:2
    for iData=1:length(NData)
        % TempData=zscore(double(NData{iData}(iscell,:)),0,2);
        TempData=double(NData{iData}(iscell(iCell),:));
        TempData=AmpNormalize(TempData,[0 100]);

            for iLaser=1:length(LaserG)
            % for iLaser=1:8

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
                   subplotLU(1,length(NData),1,iData,P);hold on
                   BaseLine=repmat(mean(tempPSTH(1:PreImgN,:),1),size(tempPSTH,1),1);
                   tempPSTH=tempPSTH-BaseLine;


                   error_area(1:size(tempPSTH,1),mean(tempPSTH,2),ste(tempPSTH')',colorLaser(iLaser,:),0.1);
                   % if iRep==1&&iLaser==1
                   %    text(15,max(mean(tempPSTH,2)),Nlabel{iData})
                   %    % set(gca,'ylim',[0 0.5])
                   % end
                   set(gca,'ylim',[-0.02 0.12])

                   % plot(PSTHLaser(:,iLaser),'color',colorLaser(iLaser,:));
                   % set(gca,'xtick',[0 PreImgN+0.5 PreImgN+0.5+PostImgN],'xticklabel',{'Start' num2str(PreImgN) num2str(PreImgN+PostImgN)})
                   % set(gca,'tickdir','out','clim',ClimScale(iData,:))
                   if iRep==3
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
     colormap(colorLaser)    
     b = colorbar;
     set(b,'position',[0.92 0.2 0.01 0.5]);
     ylabel(b, 'PV power');
     set(b,'ytick',[1:length(LaserG)]/length(LaserG),'yticklabel',LaserG);

    end
     papersizePX=[0 0 length(NData)*6 5 ];
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







