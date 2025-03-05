clear all
ProcessFolder='F:\LuSLMOnlineTest\MouseMei03\05072024\1X\';
load('C:\Users\zhangl33\Projects\GenMatCode\Plotfun\Color\colorMapPN3.mat');

ResultFolder=[ProcessFolder 'Result\'];
mkdir(ResultFolder)

TiffFolder=dir([ProcessFolder 'TSeries*'])
for i=1:size(TiffFolder,1)
    Session(i)=str2num(TiffFolder(i,1).name(end-2:end));
end

SessInfo=readtable('F:\LuSLMOnlineTest\MouseMei03\05072024\Session_Data1X.xlsx')
[~,I1,I2]=intersect(Session,SessInfo.Session);
SessInfoNeed=SessInfo(I2,:);


LastFrame=cumsum(SessInfoNeed.Repetitions.*SessInfoNeed.RepeatTimes);
FirstFrame=[1;LastFrame(1:end-1)+1];



Fall=load('F:\LuSLMOnlineTest\MouseMei03\05072024\1X\OldParam\suite2p\combined\Fall.mat');
cellInfo=Suite2pCellInfo(Fall)

PreImgN=10;
PostImgN=40;
clear IndStart IndEnd
for i=1:length(FirstFrame)
    IndStart(i,1)=FirstFrame(i)-PreImgN;
    IndEnd(i,1)=FirstFrame(i)+PostImgN-1;
end
SessInfoNeed.PreFrame=IndStart;
SessInfoNeed.PostFrame=IndEnd;
% SessInfoNeed(SessInfoNeed.PreFrame<0,:)=[];
NonSLMInd=find(SessInfoNeed.SLM==0);
SLMInd=find(SessInfoNeed.SLM==1);

SingleRep=find(SessInfoNeed.RepeatTimes==1);
MulRep=find(SessInfoNeed.RepeatTimes>1);



LaserG=unique(SessInfoNeed.LaserPower)
RepG=unique(SessInfoNeed.Repetitions)
% PmtG=unique(SessInfoNeed.PMTLevel)

deltaFoF=F2deltaFoF(Fall.F,Fall.Fneu,Fall.ops.fs);


NData={deltaFoF', Fall.spks};
Nlabel={'DeltaF', 'Spks'}

cellInfo=Suite2pCellInfo(Fall)

planeG=unique(cellInfo.iplane);
for iplane=1:length(planeG)
    PlaneC(iplane)=max(find(cellInfo.iplane==planeG(iplane)));
end


iscell=find(Fall.iscell(:,1)==1);
ClimScale=[-5 5;-20 20]
ClimScale=[-3 3;-0.2 0.2]
ClimScale=[-5 5;-10 10]

   P.xLeft=0.03;        %%%%%%Left Margin
   P.xRight=0.04;       %%%%%%Right Margin
   P.yTop=0.02;         %%%%%%Top Margin
   P.yBottom=0.1;      %%%%%%Bottom Margin
   P.xInt=0.01;         %%%%%%Width-interval between subplots
   P.yInt=0.03;         %%%%%%Height-interval between subplots



    % I1=intersect(NonSLMInd,SingleRep);
    % figure;
    % for iData=1:length(NData)
    %     % TempData=zscore(double(NData{iData}(iscell,:)),0,2);
    %     TempData=double(NData{iData}(iscell,:));
    %     for iCell=1:size(TempData)
    %         TempData(iCell,:)=AmpNormalize(TempData(iCell,:),[0 100]);
    %     end
    %     % 
    % 
    %     % for iCell=1:size(cellInfo,1)
    %         for iLaser=1:length(LaserG)
    % 
    %             I2=find(SessInfoNeed.LaserPower==LaserG(iLaser));
    %             I3=intersect(I1,I2);
    %             if ~isempty(I3)
    %                tempPSTH=[];
    % 
    %                for iSess = 2:length(I3)
    %                    TempI=SessInfoNeed.PreFrame(I3(iSess)):SessInfoNeed.PostFrame(I3(iSess));
    %                    tempPSTH(:,:,iSess)=TempData(:,TempI);
    %                end
    %                tempPSTH=squeeze(mean(tempPSTH,3));
    %                    BaseLine=repmat(mean(tempPSTH(1:PreImgN,:),1),size(tempPSTH,1),1);
    %                tempPSTH=tempPSTH-BaseLine;
    %                subplotLU(length(NData),length(LaserG),iData,iLaser,P);
    %                imagesc(tempPSTH);hold on;
    %                if iData==2
    %                set(gca,'xlim',[0 PreImgN+0.5+PostImgN] ,'xtick',[0 PreImgN+0.5 PreImgN+0.5+PostImgN],'xticklabel',{['-' num2str(PreImgN)] '0' num2str(PostImgN)})
    %                xlabel(['PV-Power ' num2str(LaserG(iLaser))])
    %                else
    %                set(gca,'xlim',[0 PreImgN+0.5+PostImgN] ,'xtick',[0 PreImgN+0.5 PreImgN+0.5+PostImgN],'xticklabel',[])
    % 
    %                end
    %                set(gca,'tickdir','out','clim',ClimScale(iData,:))
    %                colormap(ColorPN3)
    %                if iLaser >1 
    %                   set(gca,'ytick',[]);
    %                elseif   iLaser == 1 
    %                   ylabel('CellsID');
    %                   set(gca,'ytick',PlaneC+0.5,'yticklabel',{'Plane1' 'Plane2' 'Plane3'})
    % 
    %                else
    %                end
    % 
    %                for iplane=1:length(PlaneC)
    %                    plot([0 PreImgN+0.5+PostImgN],zeros(1,2)+PlaneC(iplane)+0.5,'g-')
    %                end
    %                plot([PreImgN PreImgN]+0.5,[0 length(iscell)],'k-')
    % 
    %             end
    % 
    %         end
    %      b = colorbar;
    %      set(b,'position',[0.97 0.7-(iData-1)*0.5 0.002 0.2]);
    %      ylabel(b, Nlabel{iData});
    % 
    %     % end
    % end
    %  papersizePX=[0 0 length(LaserG)*5 length(NData)*5];
    %  set(gcf, 'PaperUnits', 'centimeters');
    %  set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
    %  saveas(gcf,[ResultFolder 'AllNeuroSessReptition' num2str(RepG(iRep))],'png');
    % close all

LaserGall=[159 180 200 221 239]

PlaneNumStart=[1 PlaneC(1:end-1)+1]
jet=colormap('jet');
colorLaser=jet(1:size(jet,1)/length(LaserGall):size(jet,1),:);
colorLaser=colorLaser(3:5,:);
close all

I1=intersect(NonSLMInd,SingleRep);
figure;
for iData=1:length(NData)
        % TempData=zscore(double(NData{iData}(iscell,:)),0,2);
        for iPlane = 1:length(PlaneC)

        TempData=double(NData{iData}(iscell(PlaneNumStart(iPlane):PlaneC(iPlane)),:));
        for iCell=1:size(TempData)
            TempData(iCell,:)=AmpNormalize(TempData(iCell,:),[0 100]);
        end
        % 
        subplotLU(length(NData),length(PlaneC),iData,iPlane,P);

        % for iCell=1:size(cellInfo,1)
            for iLaser=1:length(LaserG)
               
                I2=find(SessInfoNeed.LaserPower==LaserG(iLaser));
                I3=intersect(I1,I2);
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
                   error_area(1:size(tempPSTH,2),mean(tempPSTH,1),ste(tempPSTH),colorLaser(iLaser,:),0.5);

                   % imagesc(tempPSTH);hold on;
                   % if iData==2
                   % set(gca,'xlim',[0 PreImgN+0.5+PostImgN] ,'xtick',[0 PreImgN+0.5 PreImgN+0.5+PostImgN],'xticklabel',{['-' num2str(PreImgN)] '0' num2str(PostImgN)})
                   % xlabel(['PV-Power ' num2str(LaserG(iLaser))])
                   % else
                   % set(gca,'xlim',[0 PreImgN+0.5+PostImgN] ,'xtick',[0 PreImgN+0.5 PreImgN+0.5+PostImgN],'xticklabel',[])
                   % 
                   % end
                   % set(gca,'tickdir','out','clim',ClimScale(iData,:))
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
                   plot([PreImgN PreImgN]+0.5,[0 0.05],'k-')

                end

            end
        end
         % colormap(colorLaser)
     colormap(colorLaser)    
     b = colorbar;
     set(b,'position',[0.9 0.2 0.03 0.5]);
     ylabel(b, 'PV power');
     set(b,'ytick',[1:length(LaserG)]/length(LaserG),'yticklabel',LaserG);

        % end
    end




I1=intersect(NonSLMInd,MulRep);
% close all
figure;
for iData=1:length(NData)
        % TempData=zscore(double(NData{iData}(iscell,:)),0,2);
        for iPlane = 1:length(PlaneC)

        TempData=double(NData{iData}(iscell(PlaneNumStart(iPlane):PlaneC(iPlane)),:));
        for iCell=1:size(TempData)
            TempData(iCell,:)=AmpNormalize(TempData(iCell,:),[0 100]);
        end
        % 
        subplotLU(length(NData),length(PlaneC),iData,iPlane,P);

        % for iCell=1:size(cellInfo,1)
            for iLaser=1:length(LaserG)
               
                I2=find(SessInfoNeed.LaserPower==LaserG(iLaser));
                I3=intersect(I1,I2);
                if ~isempty(I3)
                   tempPSTH=[];
                   iCount=0;
                   for iSess = 1:length(I3)
                       TempI=SessInfoNeed.PreFrame(I3(iSess)):SessInfoNeed.PostFrame(I3(iSess))
                       SessStart=SessInfoNeed.PreFrame(I3(iSess))+PreImgN;
                       for iRep = 2:SessInfoNeed.RepeatTimes(I3(iSess))
                           iCount=iCount+1;
                           S1=SessStart+(iRep-1)*SessInfoNeed.Repetitions;
                           TempI=S1-PreImgN:S1+PostImgN-1;
                           tempPSTH(:,:,iCount)=TempData(:,TempI);
                       end
                   end

                   % for iCell=1:size(tempPSTH,1)
                   %     tempPSTH(iCell,:,:)=AmpNormalize(tempPSTH(iCell,:,:),[0 100]);
                   % end
                   tempPSTH=squeeze(mean(tempPSTH,3));
                   error_area(1:size(tempPSTH,2),mean(tempPSTH,1),ste(tempPSTH),colorLaser(iLaser,:),0.5);

                   % imagesc(tempPSTH);hold on;
                   % if iData==2
                   % set(gca,'xlim',[0 PreImgN+0.5+PostImgN] ,'xtick',[0 PreImgN+0.5 PreImgN+0.5+PostImgN],'xticklabel',{['-' num2str(PreImgN)] '0' num2str(PostImgN)})
                   % xlabel(['PV-Power ' num2str(LaserG(iLaser))])
                   % else
                   % set(gca,'xlim',[0 PreImgN+0.5+PostImgN] ,'xtick',[0 PreImgN+0.5 PreImgN+0.5+PostImgN],'xticklabel',[])
                   % 
                   % end
                   % set(gca,'tickdir','out','clim',ClimScale(iData,:))
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
                   plot([PreImgN PreImgN]+0.5,[0 0.05],'k-')

                end

            end
        end
         % colormap(colorLaser)
     colormap(colorLaser)    
     b = colorbar;
     set(b,'position',[0.9 0.2 0.03 0.5]);
     ylabel(b, 'PV power');
     set(b,'ytick',[1:length(LaserG)]/length(LaserG),'yticklabel',LaserG);

        % end
    end



% NData={Fall.F, deltaFoF', Fall.spks};
% Nlabel={'F','DeltaF', 'Spks'}
% 



for iData=1:length(NData)
figure;
hold on;
    for iCell=1:size(cellInfo,1)
        % TempData=zscore(double(NData{iData}(iscell,:)),0,2);
        TempData=double(NData{iData}(iscell(iCell),:));
        TempData=AmpNormalize(TempData,[0 100]);
        plot(TempData-iCell*1,'color',[0.2 0.2 0.2]);
        EndPlane=ismember(iCell,PlaneC);
        if EndPlane
           plot([0 length(TempData)],zeros(1,2)-iCell,'g-');
        end
    end
    for iS=1:size(SessInfoNeed,1)
       [~,Laserj]=ismember(SessInfoNeed.LaserPower(iS),LaserG);

        SessStart=SessInfoNeed.PreFrame(iS)+PreImgN;
         for iRep = 1:SessInfoNeed.RepeatTimes(iS)
             S1=SessStart+(iRep-1)*SessInfoNeed.Repetitions(iS);
             plot(zeros(1,2)+S1+PreImgN,[-iCell 0],'color',colorLaser(Laserj,:));
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
     % close all
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
    I1=intersect(NonSLMInd,MulRep);

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
                   iCount=0;
                   for iSess = 1:length(I3)
                       TempI=SessInfoNeed.PreFrame(I3(iSess)):SessInfoNeed.PostFrame(I3(iSess))
                       SessStart=SessInfoNeed.PreFrame(I3(iSess))+PreImgN;
                       for iRep = 2:SessInfoNeed.RepeatTimes(I3(iSess))
                           iCount=iCount+1;
                           S1=SessStart+(iRep-1)*SessInfoNeed.Repetitions(I3(iSess));
                           TempI=S1-PreImgN:S1+PostImgN-1;
                           tempPSTH(:,iCount)=TempData(:,TempI);
                       end
                   end


                   % PSTHLaser(:,iLaser)=squeeze(mean(tempPSTH,2));
                   % error_area(1:size(tempPSTH,1),mean(tempPSTH,2),std(tempPSTH,0,2),colorLaser(iLaser,:),0.5)
                   subplotLU(1,length(NData),1,iData,P);hold on
                   BaseLine=repmat(mean(tempPSTH(1:PreImgN,:),1),size(tempPSTH,1),1);
                   tempPSTH=tempPSTH-BaseLine;

                   hold on;
                   error_area(1:size(tempPSTH,1),mean(tempPSTH,2),ste(tempPSTH')',colorLaser(iLaser,:),0.5);
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


                   % colormap(ColorPN3)

                end

            end
     colormap(colorLaser)    
     b = colorbar;
     set(b,'position',[0.92 0.2 0.01 0.5]);
     ylabel(b, 'PV power');
     set(b,'ytick',[1:length(LaserG)]/length(LaserG),'yticklabel',LaserG);

    end
% end
     papersizePX=[0 0 length(NData)*6 1*5 ];
     set(gcf, 'PaperUnits', 'centimeters');
     set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
     saveas(gcf,[ResultFolderCell 'Cell' num2str(iCell)],'png');
    close all


end







