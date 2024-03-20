clear all
TestFile='TSeries-03202024-0940-017'
DataFolder='E:\LuRecording\03202024\'

nPlane=3

close all

IndNeed=[13 14 15 16];
for j=1:length(IndNeed)
    Ind=IndNeed(j);
FrameTotal=27*3;

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
% FrameTotal=floor(length(BinDataAll)/Ly/Lx);
fclose(fileID);

XBin=reshape(BinDataAll,Ly,Lx,FrameTotal);

% clear Y
Y=zeros(Lx,Ly,floor(FrameTotal/nPlane),nPlane);
for i=1:nPlane
    Y(:,:,1:FrameTotal/nPlane,i)=XBin(1:Lx,1:Ly,i:nPlane:size(XBin,3));
end

MeanBin=squeeze(mean(Y(:,:,Ind,:),3));

MeanBin=permute(MeanBin,[2,1,3]);

figure
for iPlane=1:3
    subplotLU(2,3,1,iPlane)
    imagesc(SmoothDecDim(MeanTif(:,:,iPlane),1,3))
    caxis([0 1000])
    colormap("jet");
  
end




for iPlane=1:3
    subplotLU(2,3,2,iPlane)
    imagesc(SmoothDecDim(MeanBin(:,:,iPlane),1,3))
    caxis([0 1000]);
    colormap("jet");
end


% sum(sum(sum(abs(MeanTif-MeanBin))))
end

sum(sum(sum(abs(MeanTif(:,:,3)-MeanBin(:,:,1)))))