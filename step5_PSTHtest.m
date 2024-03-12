
load('C:\Users\zhangl33\Projects\GenMatCode\Plotfun\Color\colorMapPN3.mat');




DataFolder='E:\LuRecording\03082024_MouseMei03_VGlutAi203SLM\'
% SingP=[DataFolder 'SingleP\GPL.gpl']
SinglePSTHFolder=[DataFolder 'SinglePSTH\']
mkdir(SinglePSTHFolder);

SingPZ=[0 0 0 50 50 50 50 100 100]

BinFile=dir([TiffFolder '*TSeries*Laser*Point*.bin'])

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

PreInd=2:15
PostInd=17:30;
PlaneZ=[0 50 100];
MeanImgClim=[-600 600]

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
         subplot(1,2,iLaser)
         MultiMatrix3DPlotZ(PSTH{iLaser,iPoint},PlaneZ,0.9);
         caxis(MeanImgClim);
         Radius=10;
         colormap(ColorPN3);
         set(gca,'xlim',[0 512],'ylim',[0 512],'zlim',PlaneZ([1 end]),'View',[64 24],'zDir','reverse');
         plotCellCenter3D(SinglePxyz(:,:,iPoint), Radius, [0 1 0],1.5);

     
     end
      set(gcf, 'PaperUnits', 'centimeters');
      set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
      saveas(gcf,[SinglePSTHFolder 'Point' num2str(Point(iPoint))],'png'); 
      saveas(gcf,[SinglePSTHFolder 'Point' num2str(Point(iPoint))],'fig'); 

     close all
end




imagesc(Temp3(:,:,1))


SinglePxyz=[Y X SingPZ(iPoint)]
figure;
subplot('position',[0.01 0.01 0.88 0.95])
MultiMatrix3DPlotZ(Temp3,PlaneZ,0.9);
caxis(MeanImgClim);
Radius=8;
colormap(ColorPN3);
% set(gca,'xlim',[0 512],'ylim',[0 512],'zlim',[-80 80],'clim',ClimImg*10,'yDir','reverse');
set(gca,'xlim',[0 512],'ylim',[0 512],'zlim',PlaneZ([1 end]),'View',[64 24],'zDir','reverse');
plotCellCenter3D(SinglePxyz, Radius, [0 1 0],2.5);





nPlane=3;
for iPlane=1:nPlane

for iframe=1:length(TiffFile)
    X(:,:,iframe,iPlane)=imread([TiffFile(iframe).folder '\' TiffFile(iframe).name],iPlane);
end
    MeanImage(:,:,iPlane) = mean(X(:,:,:,iPlane),3);
end



BinFile=dir(.....bin);

fileID=fopen(BinFile(1).folder BinFile(1).name)

clear Y
BinData=fread(fileID,'unit16')
Ly=BinData(1);
Lx=BinData(2);
BinDataAll=BinData(3:end)
fclose(fileID);
XBin=reshape(BinDataAll,Ly,Lx,FrameTotal)

clear Y

for i=1:nPlane
    Y(:,:,1:FrameTotal/nPlane,i)=Xbin(Lx:-1:1,1:nPlane:size(XBin,3));;
end

figure

for iPlane=1:3
    subplot(1,3,iPlane)
    imagesc(MeanImage(:,:,iPlane))
    caxis([0 1000])
    colormap("gray");
  
end

figure
for iPlane=1:3
    subplot(1,3,iPlane)
    imagesc(squeeze(mean(Y(:,:,:,iPlane),3))')
    caxis([0 1000]);
    colormap(gray);
end