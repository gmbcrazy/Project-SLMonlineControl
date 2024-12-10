PrairieLink_RawDataStreamReg




ProcessFolder='E:\LuSLMOnlineTest\NoAnimalTest\12102024\'
FileID='TSeries-12102024-1153-007';

TiffFolder=[ProcessFolder FileID '\'];
BinFile=[ProcessFolder FileID '_onlineREG.bin'];
MotionFile=[ProcessFolder FileID '_ShiftsAndCorr.bin'];


output = ReadShiftsAndCorrFile(MotionFile);


nPlanes=3;
FrameID=1:100;
Data = Suite2pSingleChBin2Frame(BinFile,512, 512, nPlanes, FrameID);

TifData=FrameIndMultiTiffs2(TiffFolder,nPlanes,FrameID);


Data2=squeeze(Data(:,:,:,2));
TifData2=squeeze(TifData(:,:,:,2));

figure;
for i=1:100
    subplot(1,2,1)
    imagesc(Data2(:,:,i)');
    set(gca,'clim',[0 2500])
    title(i);

    subplot(1,2,2)
    imagesc(TifData2(:,:,i));
    set(gca,'clim',[0 2500])
    title(i);

    pause(0.1);


end

figure;
for i=1:100

    imshowpair(Data2(:,:,i)',TifData2(:,:,i),'falsecolor');
    title(i);

    pause(0.1);


end



TimeI=108:110;
MeanImg=mean(Data2(:,:,TimeI),3);
MeanImgTiff=mean(TifData2(:,:,TimeI),3);

figure;
imshowpair(MeanImg', MeanImgTiff, 'falsecolor');



figure;
refData=squeeze(TifData2(:,:,1));
for i=1:100
    subplot(1,2,1)
    imshowpair(SmoothDec(refData,1),SmoothDec(squeeze(TifData2(:,:,i)),1),'false')
%     imshowpair(refData,squeeze(TifData2(:,:,i)),'falsecolor')
    title(i);
%     pause(0.1);

    subplot(1,2,2)
    imshowpair(SmoothDec(refData,1),SmoothDec(squeeze(Data2(:,:,i)),1)','false')
%     imshowpair(refData,squeeze(TifData2(:,:,i)),'falsecolor')
    title(i);
    pause(0.2);

end



    imshowpair(refData,squeeze(TifData2(:,:,i)),'diff')



figure;
    subplot(1,2,1)
    imagesc(MeanImg');
    set(gca,'clim',[0 2500])

    subplot(1,2,2)
    imagesc(MeanImgTiff);
    set(gca,'clim',[0 2500])

MultiMatrix3DHeatmap(MeanImg)
