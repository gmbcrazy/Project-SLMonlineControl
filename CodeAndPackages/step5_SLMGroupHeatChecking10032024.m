% function [HeatMap, HeatMapCellNeighbour]=SLMtargetHeatMapAnalysis(ProcessFolder,Param,GlobalSavePath,SessinInfo,TargetZ)

clear all
PreMarkPointRepetition=25;    %<----------------------------------------------------------------------------------Edit,Frame # before SLM in PV
PostMarkPointRepetition=10;   %<----------------------------------------------------------------------------------Edit,Frame # after SLM in PV
PreSLMCal=15;                 %<----------------------------------------------------------------------------------Edit,Frame # before SLM to calculate baseline map
PostSLMCal=3;                 %<----------------------------------------------------------------------------------Edit,Frame # before SLM to calculate responsive map
ConfigFolder='C:\Users\zhangl33\Projects\Project-SLMonlineControl\config\';
WorkFolder='E:\LuSLMOnlineTest\SL0838-Ai203\01292025\';
ConfigFile='CurrentSLMsetting.yml';%<----------------------------------------------------------------------------------Edit, configuration file
[~,~,~,CaData,CaDataPlane,stat,yaml,confSet]=ROIToXYZ(WorkFolder,ConfigFile);
umPerPixel=mean([yaml.umPerlPixelX yaml.umPerlPixelY]);
load('C:\Users\zhangl33\Projects\Project-SLMonlineControl\CodeAndPackages\subfun\Color\colorMapPN3.mat','colorMapPN1')

nPlane=length(yaml.Zdepth_ETL);
frameRepetion=PreMarkPointRepetition+PostMarkPointRepetition; %%Total repepitions of Z series in T series setting.


PVparam.maxFrame=nPlane*frameRepetion;
PVparam.BreakPointFrame=PreMarkPointRepetition*nPlane;
PVparam.maxFrame=nPlane*frameRepetion;

%param for calculate the PSTH heatmap for online analysis
PSTHparam.PreInd=PreMarkPointRepetition-PreSLMCal:PreMarkPointRepetition;
PSTHparam.PostInd=PreMarkPointRepetition+1:PreMarkPointRepetition+PostSLMCal;
PSTHparam.Plot=1;
PSTHparam.SmoothSD=1;
PSTHparam.ColorMap=colorMapPN1;
PSTHparam.Clim=[-400 400];

%param for xml files
% XMLparam.ProcessFolder=WorkFolder;
% XMLparam.TotalRounds=confSet.NumTrial;

ProcessFolder = Get_ExpDataFolder(WorkFolder, 'SpeedStimEdgeExc', {'Data','AllIncluded','DataSum','.gpl','.xml'})



% 

ResultFolder=SumDataFolder;


PreImgN=0; %%Frame num. before imaging start to show the previous imaging data.
PreSLMFrameN=15; %Frame num. before SLM to calculate baseline level
PostSLMFrameN=3; %Frame num. after SLM to calculate response level
TrialRepTh=4; %%Minimal num. of trials to quantify SLM response
DistTh=8;%%Maximal distance by pixels considered between target cell and SLM target center.
PixFromMedCenter=8;


[AnimalInfo, DateInfo] = extractAnimalIDandDate(ProcessFolder);
SurgeryInfo = readtable('\\nimhlabstore1.nimh.nih.gov\UFNC\FNC2\Zhang\Mouse\SurgeryTable.xlsx');

[~,ic]=ismember(AnimalInfo,SurgeryInfo.MouseID)
if ~isempty(ic)
   Genotype=SurgeryInfo.Genotyping(ic);
   Virus=SurgeryInfo.Virus(ic);
else
   Genotype={'NA'};
   Virus={'NA'};
end

if exist([ProcessFolder 'Data\'])
   DataFolder=[ProcessFolder 'Data\'];
else
   DataFolder=ProcessFolder;   
end
% MP=load([ProcessFolder 'SLMIncludedIndFromIscell.mat']);
% SingP=[DataFolder 'SingleP\GPL.gpl']
SinglePSTHFolder=[DataFolder 'SinglePSTH\']
ResultFolder=SinglePSTHFolder;
mkdir(SinglePSTHFolder);

% SingPZ=[0 0 50 50 50 100 100]
BinFile=dir([DataFolder '*TSeries*GPoint*.bin']);
BinTable = PointLaser_files(BinFile);
BinTable=removevars(BinTable,{'Laser'});

TiffTable = ExpInfoMultiTiffFolder(DataFolder);
% if ~isempty(SessinInfo)
% SessInfo=readtable([ProcessFolder 'SessionInfo.xlsx']);
% SessInfo=removevars(SessInfo,{'TotalRepetition'});
% SessInfo=removevars(SessInfo,{'Power'});
% SessInfo.Properties.VariableNames{'Session'} = 'FileID';
% % SessInfo.Properties.VariableNames{'Power'} = 'Laser';
% % SessInfo.Properties.VariableNames{'TotalRepetition'} = 'totalRepetitions';
% % SessInfo.Properties.VariableNames{'Session'} = 'FileID';
% DataList = outerjoin(TiffTable, SessInfo, 'Keys', 'FileID', 'MergeKeys', true);
% DataList(isnan(DataList.totalRepetitions),:)=[];
% elseif exist('BinTable')
% 
% DataList = outerjoin(TiffTable, BinTable, 'Keys', 'FileID', 'MergeKeys', true);
% else
% 
% 
% 
% end
FunGroupI=[];
FunGroupStimKey='Group';
for i=1:size(TiffTable,1)
    if ~isempty(TiffTable.StimID{i})
       temp=findstr(TiffTable.StimID{i}{1},FunGroupStimKey);
       if ~isempty(temp)
          FunGroupI=[FunGroupI;i];
       end
    end
end

DataList=TiffTable;
DataList=DataList(FunGroupI,:);
for i=1:size(DataList,1)
    DataList.StimID{i}=DataList.StimID{i}{1};
    DataList.StimGPLInd{i}=DataList.StimGPLInd{i}{1};
end


Param.PostImgN=median(DataList.totalRepetitions);


[NeuronPos3D,NeuronPos3DRaw,CaData,CaDataPlane,stat,yaml]=Suite2pMultiPlaneROIToXYZ(DataFolder,DataList.FileID(10));
[cellIDMap,CellPixCount,MedCenter,cellBoundary]=Suite2pCellIDMapFromStat(CaData.statCell,[yaml.FOVsize_PX yaml.FOVsize_PY]);
% [cellIDMap,CellPixCount,MedCenter,cellBoundary]=Suite2pCellIDMapFromStat(stat,FovSize)
% cellBoundary = CellIDMap2Boundary(cellIDMap);
yaml.Zdepth_ETL=unique(round(yaml.Zdepth_ETL));
%% In this recording file, I forgot take single image to identify the GPL Position. 
resultPaths = findAllFiles([ProcessFolder], 'GPLFunGroup.gpl');

if exist([ProcessFolder 'SLMIncludedIndFromIscell.mat'])
    MPtemp=load([ProcessFolder 'SLMIncludedIndFromIscell.mat']);
    if  isfield(MPtemp,'Pos3DNeed')
       MP3D=MPtemp.Pos3DNeed;
   else
       MP3D=MPtemp.Pos3Dneed;
   end
    [~,LocB]=ismember(ceil(MP3D(:,3)),ceil(yaml.scan_Z(1)+yaml.Zdepth_ETL));
    MP3D(:,3)=yaml.Zdepth_ETL(LocB);
    for i=1:size(MP3D,1)
         MPName{i}=['Point ' num2str(i)];
    end
elseif ~isempty(resultPaths)
    gplFile=resultPaths{1};
    tbl=gpl2Table(gplFile);
    % XYPos=gplXYtoPixel(tbl,yaml);
    [MPPos,Z]=gplXYtoPixel(tbl,yaml);
    MPName=tbl.Name;
    MP3D=[MPPos yaml.Zdepth_ETL(Z(:,2))'];
else
    disp('No MarkPoint is available')
end

%% Load Suite2p
% suite2pPath = findAllFolders(DataFolder, 'suite2p');
% if length(suite2pPath)==1
%    CombinedSuite2p= findAllFiles([suite2pPath{1} 'combined\'], 'Fall.mat');
% end
% if length(CombinedSuite2p)==1
% Fall=load(CombinedSuite2p{1});
% cellInfo=Suite2pCellInfo(Fall);
% end
% 
% 
% if sum(DataList.totalRepetitions)~=size(Fall.F,2)
%    disp('Warning!!! Size of Frames from Suite2p does NOT match tiff nums in Tiff Folders');
% end
% SinglePxyz=[];
% SinglePxyz=Pos3DNeed;
% convertTableEntries(DataList)
% numericTable = convertStringsToNumbers(DataList)
% tableWithNumeric = convertNumericStrings(DataList);

DataList.Laser=str2num(DataList.Laser);

Laser=unique(DataList.Laser);
[StimGroup,GroupTempI]=unique(DataList.StimID);
Laser(isnan(Laser))=[];
StimGroup(isnan(StimGroup))=[];
% % SinglePxyz=Pos3DNeed(Point,:);
% SLMframe=median(DataList.PostSLMframe)-1;
SLMframe=median(DataList.PostSLMframe);


for iGroup=1:length(StimGroup)
    tempI1=DataList.StimGPLInd{GroupTempI(iGroup)};
    SLM3DGroup{iGroup}=MP3D(tempI1,:);
    SLMName{iGroup}=StimGroup{iGroup};
end

% DataList.PostSLMframe(1);

% TrialRepTh=4;
% PreImgN=10;
% PostImgN=min(90,median(DataList.totalRepetitions));
PostImgN=median(DataList.totalRepetitions);

% PostSLMFrameN=3;
% PreSLMFrameN=20;


SessInfoNeed=DataList;
LastFrame=cumsum([SessInfoNeed.totalRepetitions]);
FirstFrame=[1;LastFrame(1:end-1)+1];
clear IndStart IndEnd
for i=1:length(FirstFrame)
    IndStart(i,1)=FirstFrame(i)-PreImgN;
    IndEnd(i,1)=FirstFrame(i)+PostImgN-1;
end

SessInfoNeed.FirstFrame=FirstFrame;
SessInfoNeed.PreFrame=IndStart;
SessInfoNeed.PostFrame=IndEnd;
SessInfoNeed(SessInfoNeed.PreFrame<0,:)=[];
% NonSLMInd=find(SessInfoNeed.Laser==0&SessInfoNeed.Point>0);
% SLMInd=find(SessInfoNeed.Laser(:,1)>1&SessInfoNeed.Point>0);     %% Point  = 0 refers no SLM, < 0 refers to Group Stimuli. 
% MarkPoint=unique(SessInfoNeed.Point(SLMInd));




%% Plot SLM target location, Cell ROIs with mean imaging. 
colorGroup=jet(length(StimGroup));
colorGroup = colorGroup(randperm(length(StimGroup)),:);
MeanImgClim=[0 1];




% CellID=1:size(NeuronPos3D,1);


%% Fine Cell ROI close to SLM target, as responsive cell 
% [SLMtarget,SLMtargetCellDist]=SLMtargetMatchCell(SLM3D,NeuronPos3D,DistTh);
LaserG=Laser;
% deltaFoF=F2deltaFoF(Fall.F,Fall.Fneu,Fall.ops.fs);

colorLaser=colormap("jet");
colorLaser=colorLaser(1:size(colorLaser,1)/length(LaserG):size(colorLaser,1),:);
close
Param.Clim=[-400 400]

%% Ave SLM targets illustration
close all
ResultFolderCell=[ResultFolder 'SLMtarget\'];
mkdir(ResultFolderCell);
% PostSLMFrameN=3;
% PreSLMFrameN=30;
nPlane=length(CaDataPlane);

tempSLM3D=SLM3DGroup{iGroup}
[~,SLMplane]=ismember(tempSLM3D(:,3),yaml.Zdepth_ETL);
Param.ColorMap=colorMapPN1;
% SLMTable = [];
LaserG=setdiff(LaserG,0);
for iGroup=1:length(SLM3DGroup)
    % iCell=SLMtarget(jCell);
    tempSLM3D=SLM3DGroup{iGroup};
    [~,SLMplane]=ismember(tempSLM3D(:,3),yaml.Zdepth_ETL);
    
    % if iCell<0
    %    SLMTable(end+1,:) =  [Point(jCell) 0 0 0 1 0 PreSLMFrameN PostSLMFrameN 0 0 0];
    %    continue;
    % end
    % 
     figure;


        % TempData=double(NData{iData}(iscell(iCell),:));
        % TempData=AmpNormalize(TempData,[0 100]);

     for iLaser=1:length(LaserG)
            % for iLaser=1:8
                if LaserG(iLaser)==0
                   I2=find(SessInfoNeed.Laser(:,1)==LaserG(iLaser));      
                else
                   I2=find(SessInfoNeed.Laser(:,1)==LaserG(iLaser)&ismember(SessInfoNeed.StimID,SLMName{iGroup}));
                end

                I3=I2;


                % if ~isempty(I3)
                if length(I3)>=TrialRepTh

                   tempPSTH=[];
                   nTrial=0;
                   for iSess = 1:length(I3)
                       % TempI=SessInfoNeed.PreFrame(I3(iSess)):SessInfoNeed.PostFrame(I3(iSess));
                       % tempPSTH(:,iSess)=TempData(TempI);


                       for iS=1:length(SLMframe)
                       TempI1=SLMframe(iS);    
                       % TempI1=SessInfoNeed.FirstFrame(I3(iSess))+SLMframe(iS)-1;
                       PostStim=TempI1:TempI1+PostSLMFrameN-1;
                       PreStim=TempI1-PreSLMFrameN:TempI1-1;


                       % preSLM=[preSLM;TempData(PreStim)'];
                       % postSLM=[postSLM;TempData(PostStim)'];

                       preMap=MeanFrameIndTifFolder([DataFolder SessInfoNeed.name{I3(iSess)} '\'],nPlane,PreStim);
                       postMap=MeanFrameIndTifFolder([DataFolder SessInfoNeed.name{I3(iSess)} '\'],nPlane,PostStim);
                      
                       preMap=SmoothDecDim3(preMap,[1 1]);
                       postMap=SmoothDecDim3(postMap,[1 1]);
                           if iS==1&&iSess == 1
                              tempHeatMap=postMap-preMap;
                              tempBaseMap=preMap;
                           else
                              tempHeatMap=tempHeatMap+postMap-preMap;
                              tempBaseMap=tempBaseMap+preMap;
                           end
                           nTrial=nTrial+1;
                       end

                       % Temp1(:,:,:,iTrial)=MeanFrameIndTifFolder([DataFolder SessInfoNeed.name{I3(iSess)} '\'],nPlane,PreInd);
                       % Temp2(:,:,:,iTrial)=MeanFrameIndMultiTiffs([DataFolder SessInfoNeed.name{I3(iSess)} '\'],nPlane,PostInd);
                       % 
                   end

                   % tempHeatMap=tempHeatMap./tempBaseMap;
                   % HeatMap{jCell}{iLaser}=tempHeatMap/nTrial;
                   tempHeatMap=tempHeatMap/nTrial;

                   HeatMap{iGroup}{iLaser}=permute(tempHeatMap,[2 1 3]);
                   % HeatMap{jCell}{iLaser}=SmoothDecDim3(permute(tempHeatMap/nTrial,[2 1 3]),[1 1]);

                   % subplotLU(length(LaserG),2,iLaser,1)
                   % str=['Group' num2str(iGroup) 'Laser' num2str(LaserG(iLaser))];
                   % title([num2str(jCell) ])
                   % roiNeighbour=Center2neighbour(squeeze(HeatMap{jCell}{iLaser}(:,:,SLMplane(jCell))),SLM3D(jCell,:),PixFromMedCenter);
                   % 
                   % HeatMapCellNeighbour{jCell}{iLaser}= roiNeighbour;
                   % 
                   % imagesc(roiNeighbour');
                   % set(gca,'clim',Param.Clim);
                   % % colormap(Param.ColorMap);
                   % set(gca,'xtick',[],'ytick',[]);

                   % subplotLU(length(LaserG),2,iLaser,2)
                   % MultPlaneIs2DShow1Plane(HeatMap{iGroup}{iLaser}, [], tempSLM3D, [], yaml.Zdepth_ETL, SLMplane, [0 1 0], Param.Clim);
                   MultiPlanes2DShow(HeatMap{iGroup}{iLaser}, [], tempSLM3D, [], yaml.Zdepth_ETL, [0 1 0], Param.Clim);

                   % imagesc(roiNeighbour);
                   % set(gca,'clim',Param.Clim);
                   colormap(Param.ColorMap);
                   set(gca,'xtick',[],'ytick',[]);


                end

      end
      % title(['MP' num2str(num2str(MarkPoint(jCell))) 'Cell' num2str(iCell)]);
      papersizePX=[0 0 3*8 length(LaserG)*8 ];
      set(gcf, 'PaperUnits', 'centimeters');
      set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
      saveas(gcf,[ResultFolderCell 'SLMTargetFunGroupMap' num2str(iGroup)],'png');
     close all


end





