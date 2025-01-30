

%% Load Data
clear all
% ProcessFolder='F:\LuSLMOnlineTest\04222024\SingleP\30PixelFromEdgeExc\';
load('C:\Users\User\Project-SLMonlineControl\subfun\Color\colorMapPN3.mat');
ConfigFolder='C:\Users\User\Project-SLMonlineControl\config\';
WorkingFolder='E:\LuSLMOnlineTest\SL0777-Ai203\12182024\';%<--------------------     --Edit, Data folder
PreDefTseriesFolder=[CConfigFolder 'PreGenerateTseriesMultiZ\'];



ConfigFile='SLMsetting.yml';%<---------------------------------------------------------Edit, configuration file
[~,~,~,CaData,CaDataPlane,stat,yaml,confSet]=ROIToXYZ(ConfigFolder);
umPerPixel=mean([yaml.umPerlPixelX yaml.umPerlPixelY]);
ProcessFolder=[WorkingFolder 'SingleP\Top12SpeedStimEdgeExc\'];%<----------------------Edit, Data folder
SumDataFolder=[ProcessFolder '\DataSum\'];
DataLogFolder=[ProcessFolder 'DataLog\'];

load([ProcessFolder 'SLMFunGroup.mat'],'Group','FinalPos3D','FinalCellstat','FinalFunScore','confSetFinal','SLMTableOrigin','SLMTable','ROIparam','SLMRes','sampleN','SLMTestParam','SLMIncludedIndFromIscell','FunScore','yaml','Cellstat');
load([ConfigFolder 'SpontBeh5T_Z11Frame550.mat'])




% PreMarkPointRepetition=25;    %<----------------------------------------------------------------------------------Edit,Frame # before SLM in PV
% PostMarkPointRepetition=10;   %<----------------------------------------------------------------------------------Edit,Frame # after SLM in PV
PreSLMCal=15;                   %<----------------------------------------------------------------------------------Edit,Frame # before SLM to calculate baseline map
PostSLMCal=3;                   %<----------------------------------------------------------------------------------Edit,Frame # before SLM to calculate responsive map

nPlane=length(confSet.ETL);

PVparam.InterMPRepetition=[40 40 80 80 30 50 40 60 90 50 30];
frameRepetition=sum(PVparam.InterMPRepetition); %%Total repepitions of Z series in T series setting;



PVparam.maxFrame=nPlane*frameRepetition;
PVparam.BreakPointFrame=PVparam.InterMPRepetition(1:end-1)*nPlane;
% PVparam.InterMPFrame=[40 60 30 20]*nPlane;
PVparam.TrialMPSwitch=length(PVparam.InterMPRepetition)-1;
PVparam.nPlane=nPlane;




%param for calculate the PSTH heatmap for online analysis
% PSTHparam.PreInd=PreMarkPointRepetition-PreSLMCal:PreMarkPointRepetition;
% PSTHparam.PostInd=PreMarkPointRepetition+1:PreMarkPointRepetition+PostSLMCal;
% PSTHparam.Plot=1;
PSTHparam.SmoothSD=1;
PSTHparam.ColorMap=colorMapPN1;
PSTHparam.Clim=[-400 400];

%param for xml files
% XMLparam.Laser=1.5;               %<-------------------------------------------------------------------------------Edit, starting laser power to test    
% XMLparam.RoundID=1;               %starting round

XMLparam.ProcessFolder=ProcessFolder;
XMLparam.TotalRounds=confSet.NumTrial;
% PointAll=1:size(Pos3Dneed,1);


%param for ROI neighbourhood to determine wether there is SLM response.
% ROIparam.TotalSLMPos3D=Pos3Dneed;    %%such that ROIparam.PointAll=1:size(Pos3Dneed,1)
% ROIparam.PointAll=PointAll;
% ROIparam.PlaneZ=PlaneZ;
% ROIparam.CellSize=20;                %%normal neuron diameter by um;        
% ROIparam.threshold_percentage=0.3;   %%thereshold to define responsive fields SLM responsive heatmap: percentage*Peak rate
% ROIparam.thNum=15;                   %%Minimal single responsive field by pixels
% ROIparam.max_distance=ceil(ROIparam.CellSize*2/3/umPerPixel);  %% 2/3 diameter of a cell by pixel as maximal response region-SLM center distance
% ROIparam.min_region_size=10;
% ROIparam.PeakTh=200;
% ROIparam.min_merged_region_size=20;  %%Minimal total size of responsive fields by pixels
% ROIparam.contourMethod='perim';      %%Method to detect ROI boader of responsive fields
% ROIparam.NeighbourHfWidthPixel=20;   %%PixFromMedCenter: Number of pixels from the median center of each ROI to get the ROI neighborhood.
% ROIparam.umPerPixel=mean([yaml.umPerlPixelX yaml.umPerlPixelY]);
% ROIparam.Colormap=colorMapPN1;                  
% ROIparam.LaserPower=confSet.UncagingLaserPower;   

XMLparam.ShamPossibility=0.2;
XMLparam.SwitchXMLPostMPFrame=10;
XMLparam.ProcessFolder=ProcessFolder;



% PreMarkPointRepetition=[40 60 50];
% PostMarkPointRepetition=20;

%% 
numGPUs=0;      %%Do not use GPU, assume in general the aquisition PC has no GPU. 
FileType=2;   %Choose a specific bin file as reference for motion correction
% ProcessFolder='E:\LuSLMOnlineTest\SL0777-Ai203\12122024\'
% RefFile=[DataFolder 'TSeries-12132024-1247-023.bin'];
% [RegOps, RegImg] = LoadRegRefFile(RefFile, FileType, numGPUs, [512,512,3,30]);
% 
% FileType=0;   %Choose a pre-recorded multi-tif files for motion correction
% RefFile=[];
% RefFile='E:\LuSLMOnlineTest\SL0777-Ai203\12122024\TSeries-12122024-0938-000\';
% 
% [RegOps, RegImg] = LoadRegRefFile(RefFile, FileType);
% MultiMatrix3DHeatmap(RegImg)
numGPUs=0;      %%Do not use GPU, assume in general the aquisition PC has no GPU. 
FileType=2;   %Choose a specific bin file as reference for motion correction
% ProcessFolder='E:\LuSLMOnlineTest\SL0777-Ai203\12122024\'
% RefFile=[DataFolder 'TSeries-12132024-1247-023.bin'];
% [RegOps, RegImg] = LoadRegRefFile(RefFile, FileType, numGPUs, [512,512,3,30]);
% 
FileType=0;   %Choose a pre-recorded multi-tif files for motion correction
RefFile=[];
RefFile=[WorkingFolder 'RegRef2\'];
% 
[RegOps, RegImg] = LoadRegRefFile(RefFile, FileType,numGPUs);

% FileType=1;   %Choose suite2p folder, using ops.meanImg for motion correction
% RefFile=[WorkingFolder 'suite2p\'];
% [RegOps, RegImg] = LoadRegRefFile(RefFile, FileType,numGPUs);

XMLparam.DoRegistration=1;
XMLparam.RegRefOps=RegOps;
XMLparam.RegRefImg=RegImg;  
XMLparam.RegRefSource=RefFile;  

disp('Referrence Img for Motion correction updated')


%%
XMLTable=[];
PSTHmap=[];
CountExp=1;
TotalGroupIDs=[1 2 3];   %% All possible Functional Group IDs.

XMLparam.LoadGPL=1;

[~,~]=PV_LinkExcuteXMLFunGroup(XMLparam,PVparam);




for iExp=1:length(ExpFileInfo)
XMLTable=[XMLTable;tempXMLTable{iExp}];
[indexVector, stimulusIDVector, PrePostStimuliVector] = getPSTHFrames(PVparam.InterMPRepetition, PreSLMCal, PostSLMCal);
tempPSTHmap = CalMultiPSTHBin(ExpFileInfo(iExp).binFile, confSet, indexVector, stimulusIDVector, PrePostStimuliVector);
PSTHmap=cat(3,PSTHmap,tempPSTHmap);
end

CountExp=CountExp+1;
OutPng=[SumDataFolder 'CheckFunExp' num2str(CountExp) '.png'];



[averagedPSTHmap,nSample] = AveragePSTHByGroupLaser(XMLTable, PSTHmap,TotalGroupIDs);
Zdepth=confSet.scan_Z+confSet.ETL;



figure;
for iFun = 1:length(nSample)
    if iFun<length(nSample)
       FinalFunScore(Group(iFun).Indices, 1) = iFun;
       Pos3D=FinalPos3D(Group(iFun).Indices,:);
       RealMPInd=find((isnan(FinalFunScore(Group(iFun).Indices,2)))==0);
       AddMPInd=find((isnan(FinalFunScore(Group(iFun).Indices,2)))==1);

     for iplane = 1:nPlane
        if nSample(iFun)>0
           subplotLU(length(nSample), nPlane, iFun, iplane);
            
            % Display the current plane's image
            imagesc(averagedPSTHmap(:,:,iplane,iFun)');
          
            I = find(abs(Pos3D(:,3) - Zdepth(iplane)) < 0.1);
            I1=intersect(I,RealMPInd);
            I2=intersect(I,AddMPInd);
            if ~isempty(I1)
                plotCellCenter(Pos3D(I1,[2 1]), 7, [0.1 0.9 0.1],1);    
            end

            if ~isempty(I2)
                plotCellCenter(Pos3D(I2,[2 1]), 7, [0.9 0.1 0.1],1);    
            end
        end
     end
    else
     for iplane = 1:nPlane
        if nSample(iFun)>0
           subplotLU(length(nSample), size(Img,3), iFun, iplane);
            
            % Display the current plane's image
            imagesc(averagedPSTHmap(:,:,iplane,iFun)');

        end
     end




    end


end

papersizePX=[0 0 nPlane*4 length(nSample)*4]
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
saveas(gcf,OutPng,'png');











step5_SubStep1_GetResponse;



XMLTable=[XMLTable;tempXMLTable];
[indexVector, stimulusIDVector, PrePostStimuliVector] = getPSTHFrames(PVparam.InterMPRepetition, PreSLMCal, PostSLMCal);
tempPSTHmap = CalMultiPSTHBin(ExpFileInfo(CountExp).binFile, confSet, indexVector, stimulusIDVector, prePostStimuliVector);
PSTHmap=cat(PSTHmap,tempPSTHmap,3);
CountExp=CountExp+1;



XMLTable=[groupIDs(:) laserPowers(:)];




sizeV=size(frameCal);
if length(sizeV)==4
   nPlane=sizeV(4);
else
   nPlane=1;
end
nPlane=size(frameCal,4);

stimID=unique(stimulusIDVector);
for iStim=1:length(stimID)
    NeedI=stimulusIDVector==stimID(iStim);
    tempframeID=indexVector(NeedI);
    tempPrePostID=prePostStimuliVector(NeedI);
    tempPreID=tempframeID(tempPrePostID==-1);
    tempPostID=tempframeID(tempPrePostID==1);


    if nPlane==1
    preframeCal(:,:,iStim) = squeeze(mean(frameCal(:,:,tempPreID),3));
    postframeCal(:,:,iStim) = squeeze(mean(frameCal(:,:,tempPostID),3));
    
    else
    preframeCal(:,:,iStim,:) = squeeze(mean(frameCal(:,:,tempPreID,:),3));
    postframeCal(:,:,iStim,:) = squeeze(mean(frameCal(:,:,tempPostID,:),3));
    end
end


function PSTHtemp=PSTHmapCal(ExpFileInfo(CountExp).binFile,PSTHparam,confSet)


         PreData = Suite2pSingleChBin2Frame(BinFile, confSet.SLM_Pixels_Y, confSet.SLM_Pixels_X, length(confSet.ETL), PSTHparam.PreInd);
         PostData= Suite2pSingleChBin2Frame(BinFile, confSet.SLM_Pixels_Y, confSet.SLM_Pixels_X, length(confSet.ETL), PSTHparam.PostInd);
         PSTHtemp=squeeze(mean(PostData,3)-mean(PreData,3));
         if PSTHparam.SmoothSD>0
         PSTHtemp=SmoothDecDim3(PSTHtemp,PSTHparam.SmoothSD);
         end
end



PSTHmap(:,:,:,ixml)=PSTHmapCal(ExpFileInfo(CountExp).binFile,PSTHparam,confSet);



SLMTrialMap=




[tempXMLTable,ExpFileInfo(CountExp)]=PV_LinkExcuteXMLFunGroup(XMLparam,PVparam);










