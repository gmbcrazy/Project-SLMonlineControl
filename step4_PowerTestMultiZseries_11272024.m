clear all
% TestFile='TSeries-04222024-0926-040'

load('C:\Users\User\Project-SLMonlineControl\subfun\Color\colorMapPN3.mat');
ConfigFolder='C:\Users\User\Project-SLMonlineControl\config\';
SLMsettingFile='SLMsetting.yml';
confSet = ReadYaml([ConfigFolder '\' SLMsettingFile]);

nPlane=length(confSet.ETL)
% DataFolder='F:\LuSLMOnlineTest\04222024\Data\'
ProcessFolder='E:\LuSLMOnlineTest\SL0777-Ai203\12032024\SingleP\Top5SpeedStimEdgeExc\';
DataFolder=[ProcessFolder 'Data\'];
mkdir(DataFolder);
DataLogFolder=[ProcessFolder 'DataLog\'];

Temp1=load([ProcessFolder 'SLMIncludedIndFromIscell.mat'],'Pos3Dneed');
AllTestPoints3D=Temp1.Pos3Dneed; clear Temp1
iCount=1;

%%



XMLparam.SwitchXMLPostMPFrame=6;                   %%<-----------MarkPoint switching occurs after 10 Repetitions of nplanes of Zseries.
XMLparam.ProcessFolder=ProcessFolder;

XMLparam.RoundID=6;
XMLparam.PointList=[1:10];                     %%<-----------nP, NumOfTestedPoints, nP + 1 = NumOfZseries 
XMLparam.Laser=[repmat(1.6,1,10)] ;                  %%<-----------laser values, could be 1 value of a vector with length of nP



PowerTestPVPar.nPlane=nPlane;            
PowerTestPVPar.ZRepetition=31;                      %%<-----------NumOfRepeition in each Zseries of Tseries in PV
PowerTestPVPar.Ziteration=11;                        %%<-----------NumOfZseries in Tseries in PV
PowerTestPVPar.InterMPRepetition=repmat(PowerTestPVPar.ZRepetition,1,PowerTestPVPar.Ziteration);
frameRepetition=PowerTestPVPar.ZRepetition*PowerTestPVPar.Ziteration;
PowerTestPVPar.maxFrame=nPlane*frameRepetition;
% PowerTestPVPar.BreakPointFrame=PowerTestPVPar.InterMPRepetition(1:end-1)*nPlane;
% PowerTestPVPar.InterMPFrame=[40 60 30 20]*nPlane;
% PowerTestPVPar.TrialMPSwitch=length(PowerTestPVPar.InterMPRepetition)-1;


if sum(abs([length(XMLparam.Laser) length(XMLparam.PointList)]+1-PowerTestPVPar.Ziteration))~=0
   disp('Check whether # Point, Laser levels and Zseries match')
end



[XMLTable{iCount},FileGenerateInfo(iCount)]=PV_LinkPowerTest_MultiZseries(XMLparam,PowerTestPVPar)
[checkXMLTable{iCount},UnMatchI{iCount}]=MPxmlExcuteMatchCheck(FileGenerateInfo(iCount),XMLTable{iCount},AllTestPoints3D,confSet);




iCount=iCount+1
%%



pause(30)
for iTest=1:15
    [tempXMLTable{CountExp},ExpFileInfo(CountExp)]=PV_LinkExcuteXMLFunGroup(XMLparam,PowerTestPVPar2);
    pause(3);
end


% PowerTestPVPar.maxFrame=nPlane*frameRepetition;
% PowerTestPVPar.BreakPointFrame=PowerTestPVPar.InterMPRepetition(1:end-1)*nPlane;
% % PowerTestPVPar.InterMPFrame=[40 60 30 20]*nPlane;
% PowerTestPVPar.TrialMPSwitch=length(PowerTestPVPar.InterMPRepetition)-1;
% PowerTestPVPar.nPlane=nPlane;
% 

% PV_LinkExcuteXML(XMLparam,PowerTestPVPar,confSet);

%     TestFile='TSeries-11182024-1006-042';
    TestFile=['TSeries-11182024-1006-070'];


close all
StartingFrame=210;
IndNeedTiff=[StartingFrame:StartingFrame+2];

% IndNeedTiff=[StartingFrame StartingFrame+2:StartingFrame+2+2];

IndNeedBin=[StartingFrame:StartingFrame+2]
FrameTotal=frameRepetition*3;
% FrameTotal=1650
for j=1:length(IndNeedTiff)
    Ind=IndNeedTiff(j);
    IndBin=IndNeedBin(j);

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




% DataFolder='E:\LuSLMOnlineTest\SL0777-Ai203\10302014\SingleP\Top13SpeedStimEdgeExc\Data\'
StartFile=9;
for iTest=1:10

%     TestFile='TSeries-11182024-1006-042';
    TestFile=['TSeries-11182024-1006-0' num2str(StartFile+iTest)];


close all
StartingFrame=210;
IndNeedTiff=[StartingFrame:StartingFrame+2];

% IndNeedTiff=[StartingFrame StartingFrame+2:StartingFrame+2+2];

IndNeedBin=[StartingFrame:StartingFrame+2]
FrameTotal=frameRepetition*3;
% FrameTotal=1650
for j=1:length(IndNeedTiff)
    Ind=IndNeedTiff(j);
    IndBin=IndNeedBin(j);

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
figure(2)

saveas(gcf,[DataLogFolder  TestFile '.png'],'png');
saveas(gcf,[DataLogFolder  TestFile],'fig');


end
% sum(sum(sum(abs(MeanTif(:,:,3)-MeanBin(:,:,1)))))

% 
% 40
% 
% 59
% 
% 62
% 
% 


