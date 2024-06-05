clear all
TestFile='TSeries-04222024-0926-040'
DataFolder='F:\LuSLMOnlineTest\04222024\SingleP\30PixelFromEdgeExc\Data\'

nPlane=3;
NumRepetition=80;
close all

IndNeedTiff=[60 62 63];
IndNeedBin=[60 61 62]
for j=1:length(IndNeedTiff)
    Ind=IndNeedTiff(j);
    IndBin=IndNeedBin(j);
FrameTotal=NumRepetition*3;

TiffPath=[DataFolder TestFile '\'];
BinPath=[DataFolder TestFile '*.bin'];
BinFile=dir(BinPath);
fileID=fopen([BinFile(1).folder '\' BinFile(1).name]);

ImageSeq=FrameIndMultiTiffs(TiffPath,nPlane,Ind);
MeanTif=squeeze(mean(ImageSeq,3));


clear Y
BinData=fread(fileID,'uint16');
Ly=512;
Lx=512;
% BinDataAll=BinData(:);
BinDataAll=BinData(1:Ly*Lx*FrameTotal);
% FrameTotal=floor(length(BinData)/Ly/Lx)
% BinDataAll=BinData(1:Ly*Lx*FrameTotal);

fclose(fileID);

XBin=reshape(BinDataAll,Ly,Lx,FrameTotal);

% clear Y
Y=zeros(Lx,Ly,floor(FrameTotal/nPlane),nPlane);
for i=1:nPlane
    Y(:,:,1:FrameTotal/nPlane,i)=XBin(1:Lx,1:Ly,i:nPlane:size(XBin,3));
end

MeanBin=squeeze(mean(Y(:,:,IndBin,:),3));

MeanBin=permute(MeanBin,[2,1,3]);

figure
for iPlane=1:3
%     subplotLU(2,3,1,iPlane)
    subplot(2,3,iPlane)

    imagesc(SmoothDec(MeanTif(:,:,iPlane),1))
    caxis([0 600])
    colormap("jet");
  
end




for iPlane=1:3
%     subplotLU(2,3,2,iPlane)
        subplot(2,3,iPlane+3)

    imagesc(SmoothDec(MeanBin(:,:,iPlane),1))
    caxis([0 600]);
    colormap("jet");
end


% sum(sum(sum(abs(MeanTif-MeanBin))))
end

sum(sum(sum(abs(MeanTif(:,:,3)-MeanBin(:,:,1)))))


40

59

62




