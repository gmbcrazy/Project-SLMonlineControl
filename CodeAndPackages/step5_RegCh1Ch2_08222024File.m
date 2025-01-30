
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
ProcessFolder='F:\LuSLMOnlineTest\SL0242-Ai203\08222024\';
PreImgN=10;   %%Define the frame Nums before each imaging file start to show, coming from the previous file
PostImgN=80;  %%Define the frame Nums after each imaging file start to show
ResultFolder=[ProcessFolder 'Result\'];
mkdir(ResultFolder)
SLMframe=[0 40];

load('C:\Users\zhangl33\Projects\GenMatCode\Plotfun\Color\colorMapPN3.mat');

[Pos3D,Pos3DRaw,CaData,CaDataPlane,stat,yaml]=Suite2pMultiPlaneROIToXYZ(ProcessFolder);
[cellIDMap,CellPixCount,MedCenter]=Suite2pCellIDMapFromStat(CaData.statCell,[yaml.FOVsize_PX yaml.FOVsize_PY]);
cellBoundary = CellIDMap2Boundary(cellIDMap);
yaml.Zdepth_ETL=unique(round(yaml.Zdepth_ETL));

CellID=1:size(Pos3D,1);







RedChFolder='F:\LuSLMOnlineTest\SL0242-Ai203\08222024\2ChTSeries\TSeries-08222024-0916-091\';
Fix=AmpNormalizeDim(permute(double(CaData.PlaneMeanImg),[2 1 3]),3,[1 99.9]);

for iCh=1:2
    RawImg1=FrameIndMultiTiffs2(RedChFolder,length(yaml.Zdepth_ETL),1:500,iCh);
    meanRaw{iCh}=squeeze(mean(RawImg1,3));
end
movingCh2=AmpNormalizeDim(permute(double(meanRaw{2}),[2 1 3]),3,[1 99.9]);
movingCh1=AmpNormalizeDim(permute(double(meanRaw{1}),[2 1 3]),3,[1 99.9]);


% % Perform image registration using histogram matching and demons algorithm
% Moving2 = imhistmatch(movingCh2,Fix);
% [D,MovingReg] = imregdemons(Moving2,Fix,[500 400 200],'AccumulatedFieldSmoothing',1.3);  %%D is the non-linear transform of pixel from Moving to Fixed
% 
%     %% Using the D transform to transform Moving CellID map to Fixed imaging
% MovingCellReg = imwarp(MovingCell,D,'Interp','nearest','SmoothEdges',false);
% 


clear movedCh1 movedCh2 tform
movedCh1=zeros(size(Fix));
[optimizer, metric]  = imregconfig('monomodal');
for iplane=1:3
    moving = squeeze(movingCh2(:,:,iplane));
    fix = squeeze(Fix(:,:,iplane));
tform= imregtform(moving,fix,"affine",optimizer,metric);
movedCh1(:,:,iplane) = imwarp(squeeze(movingCh1(:,:,iplane)),tform,"OutputView",imref2d(size(fix)));
movedCh2(:,:,iplane) = imwarp(squeeze(movingCh2(:,:,iplane)),tform,"OutputView",imref2d(size(fix)));
end

figure;
MultiPlanes2DShow(movedCh2, cellBoundary, Pos3D, [],yaml.Zdepth_ETL, colorCell, MeanImgClim)
papersizePX=[0 0 20 7];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
saveas(gcf,[ResultFolder 'Ch2'],'png');
close all


figure;
% MultiPlanes2DShow(movedCh1, cellBoundary, Pos3D, CellID,yaml.Zdepth_ETL, colorCell, MeanImgClim)
MeanImgClim=[0.05 0.5]
MultiPlanes2DShow(movedCh1, cellBoundary, Pos3D, [],yaml.Zdepth_ETL, colorCell, MeanImgClim)
papersizePX=[0 0 20 7];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
saveas(gcf,[ResultFolder 'Ch1'],'png');
close all



% figure;
% % ImgCh1=AmpNormalize(permute(double(meanImgCh1),[2 1 3]),[1 99]);
% MultiMatrix3DPlotZ(movedCh1,yaml.Zdepth_ETL,0.9);
% caxis(MeanImgClim);

% figure;
% 
% for iplane=1:length(yaml.Zdepth_ETL)
% ImgCh2(:,:,iplane)=CaDataPlane(iplane).ops.meanImg;
% end
% ImgCh2=AmpNormalize(permute(double(ImgCh2),[2 1 3]),[1 99]);
% 

figure;

for iplane=1:3
I=find(Pos3D(:,3)==yaml.Zdepth_ETL(iplane));
% plotCellBoundary(cellBoundary, CellID)
% imagesc(ImgCh2(:,:,iplane));hold on;
RGBImage = cat(3, movedCh1(:,:,iplane), movedCh2(:,:,iplane)*0.9, zeros(size(movedCh1(:,:,iplane))));


% Load the image
% RGBImage = imread('/mnt/data/image.png');

% Extract the red and green channels
% redChannel = AmpNormalize(log(movedCh1(:,:,iplane)+0.0000001),[0.5 99]).^4;
% greenChannel = AmpNormalize(log(movedCh2(:,:,iplane)+0.0000001),[0.5 99]);
% redChannel = SmoothDec(AmpNormalize(movedCh1(:,:,iplane),[0 100]),[1 1]);
redChannel = SmoothDec(AmpNormalize(movedCh1(:,:,iplane),[0 100]),[0.5 0.5]);

% EnRed = adapthisteq(redChannel,'NumTiles',[8 8],'ClipLimit',0.05);

% greenChannel = AmpNormalize(movedCh2(:,:,iplane),[0 100]);
greenChannel = Fix(:,:,iplane);

lowCutoff=10;
highCutoff=150;
EnRed=redChannel;
% EnRed = fftFilterImage(EnRed, lowCutoff, highCutoff);
EnRed = adapthisteq(EnRed,'NumTiles',[8 8],'ClipLimit',0.05);
% EnRed = SmoothDec(EnRed,[1 1]);

% Convert the red and green channels to binary masks using a threshold
redMask = imbinarize(EnRed, 0.75);  % Adjust threshold as needed
minPixTh=15;
redMask=MaskThreshold(redMask, minPixTh);
% imagesc(redMask)


% EnGreen = SmoothDec(fftFilterImage(adapthisteq(greenChannel), lowCutoff, highCutoff),[0.5 0.5]);
EnGreen = adapthisteq(greenChannel,'NumTiles',[8 8],'ClipLimit',0.02);
% EnGreen = SmoothDec(EnGreen,[0.5 0.5]);
% plotCellBoundary3D(cellBoundary(I), Pos3D(I,3),[0 1 0],0.5)
greenMask = imbinarize(EnGreen, 0.85);  % Adjust threshold as needed
minPixTh=40;
greenMask=MaskThreshold(greenMask, minPixTh);
% imagesc(greenMask)


overlapMask = redMask & greenMask;
minPixTh=15;
overlapMask=MaskThreshold(overlapMask, minPixTh);



overlapStats = struct2array(regionprops(overlapMask, 'Area'));
redStats = struct2array(regionprops(redMask, 'Area'));
greenStats = struct2array(regionprops(greenMask, 'Area'));


subplot(3,3,1+(iplane-1)*3);
RGBImage = cat(3, EnRed, zeros(size(EnRed)), zeros(size(EnRed)));
imshow(RGBImage);
hold on;
contour(redMask,'b')
Info=['Red Masks: ' num2str(sum(redStats)) ' Pixel with ' num2str(length(redStats)) ' Objects (Cells)'];
title(Info);

subplot(3,3,2+(iplane-1)*3);
RGBImage = cat(3, zeros(size(EnGreen)), greenChannel, zeros(size(EnGreen)));
imshow(RGBImage);hold on;
contour(greenMask,'b')
Info=['Green Masks: ' num2str(sum(greenStats)) ' Pixel with ' num2str(length(greenStats)) ' Objects (Cells)'];
title(Info);


subplot(3,3,3+(iplane-1)*3);
RGBImage = cat(3, EnRed, EnGreen, zeros(size(EnGreen)));
imshow(RGBImage);hold on;
contour(overlapMask,'b')
Info=['Green Masks: ' num2str(sum(overlapStats)) ' Pixel with ' num2str(length(overlapStats)) ' Objects (Cells)'];
title(Info);

end

papersizePX=[0 0 20 20];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
saveas(gcf,[ResultFolder 'Ch1Ch2'],'png');

