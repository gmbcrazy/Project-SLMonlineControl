clear all

load('C:\Users\zhangl33\Projects\GenMatCode\Plotfun\Color\colorMapPN3.mat');




DataFolder='E:\LuRecording\03222024Test\'
% SingP=[DataFolder 'SingleP\GPL.gpl']
SinglePSTHFolder=[DataFolder 'SinglePSTH\']
mkdir(SinglePSTHFolder);

% SingPZ=[0 0 50 50 50 100 100]

BinFile=dir([DataFolder '*TSeries*Laser*Point*.bin'])
SingPZ=[100 100]

PointFile=[];
LaserFile=[];
BinFileFolder={};
for iFile=1:length(BinFile)
    I1=findstr(BinFile(iFile).name,'Laser');
    I2=findstr(BinFile(iFile).name,'Point');
    I3=findstr(BinFile(iFile).name,'.bin');
    LaserFile(iFile)=str2num(BinFile(iFile).name(I1+5:I2-1));
    PointFile(iFile)=str2num(BinFile(iFile).name(I2+5:I3-1));
    BinFileFolder{iFile}=BinFile(iFile).name(1:I1-1);
end

Laser=unique(LaserFile)
Point=unique(PointFile)
nPlane=3
PreInd=2:25
PostInd=30:33;
PlaneZ=[0 50 100];
MeanImgClim=[-800 800]

% SingPZ=100;

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
             pointInfo=xmlInfo.PVMarkPointSeriesElements.PVMarkPointElement.PVGalvoPointElement.Point.Attributes;
             X=str2num(pointInfo.X)*512;
             Y=str2num(pointInfo.Y)*512;
             SinglePxyz(:,:,iPoint)=[Y X SingPZ(iPoint)];

             end
         end
         Temp3=mean(Temp2,4)-mean(Temp1,4);
         Temp3=SmoothDecDim(Temp3,1,3);
    
         % subplot(3,3,iPoint)
          % figure;
%          MultiMatrix3DPlotZ(Temp3,PlaneZ,0.9);
%         caxis(MeanImgClim);
%          Radius=8;
%         colormap(ColorPN3);
% % set(gca,'xlim',[0 512],'ylim',[0 512],'zlim',[-80 80],'clim',ClimImg*10,'yDir','reverse');
%         set(gca,'xlim',[0 512],'ylim',[0 512],'zlim',PlaneZ([1 end]),'View',[64 24],'zDir','reverse');
%           plotCellCenter3D(SinglePxyz, Radius, [0 1 0],2.5);

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






% SinglePxyz=[Y X SingPZ(iPoint)]
% figure;
% subplot('position',[0.01 0.01 0.88 0.95])
% MultiMatrix3DPlotZ(Temp3,PlaneZ,0.9);
% caxis(MeanImgClim);
% Radius=8;
% colormap(ColorPN3);
% % set(gca,'xlim',[0 512],'ylim',[0 512],'zlim',[-80 80],'clim',ClimImg*10,'yDir','reverse');
% set(gca,'xlim',[0 512],'ylim',[0 512],'zlim',PlaneZ([1 end]),'View',[64 24],'zDir','reverse');
% plotCellCenter3D(SinglePxyz, Radius, [0 1 0],2.5);
% 




% % nPlane=3;
% % for iPlane=1:nPlane
% % 
% % for iframe=1:length(TiffFile)
% %     X(:,:,iframe,iPlane)=imread([TiffFile(iframe).folder '\' TiffFile(iframe).name],iPlane);
% % end
% %     MeanImage(:,:,iPlane) = mean(X(:,:,:,iPlane),3);
% % end



PreInd=2:25
PostInd=30:33;

% PreInd=2:11;
% PauseInd=11;
% PostInd=13:26;
PlaneZ=[0 50 100];
clear PSTHtemp;
% BinFile=dir([DataFolder '*TSeries*Laser1.12Point6.bin'])
FrameTotal=(58-1)*nPlane
% Skip=PauseInd*3+1;




Ly=512;
Lx=512;


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

            PSTHtemp(:,:,:,iTrial)=nanmean(Y(:,:,PostInd-1,:),3)-nanmean(Y(:,:,PreInd,:),3);
             
             

         end
         PSTHtemp=squeeze(mean(PSTHtemp,4));
         PSTHtemp=permute(PSTHtemp,[2,1,3]);
         PSTHtemp=SmoothDecDim(PSTHtemp,1,3);

         PSTHBin{iLaser,iPoint}=PSTHtemp;
    end
    toc
end









close all
papersizePX=[0 0 22 18];
for iPoint = 1:length(Point)
    figure;

     for iLaser=1:length(Laser)
         subplot(length(Laser),2,1)
         MultiMatrix3DPlotZ(PSTH{iLaser,iPoint},PlaneZ,0.9);
         caxis(MeanImgClim);
         Radius=10;
         colormap(ColorPN3);
         set(gca,'xlim',[0 512],'ylim',[0 512],'zlim',PlaneZ([1 end]),'View',[64 24],'zDir','reverse');
         plotCellCenter3D(SinglePxyz(:,:,iPoint), Radius, [0 1 0],1.5);
         ylabel('Tiff Data')


         subplot(length(Laser),2,2)
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

     % close all
end





