clear all

load('C:\Users\zhangl33\Projects\GenMatCode\Plotfun\Color\colorMapPN3.mat');

% load('C:\Users\User\Project-SLMonlineControl\subfun\Color\colorMapPN3.mat');

ProcessFolder='\\nimhlabstore1.nimh.nih.gov\UFNC\UFNC_Data\DATA_LZ\Sutie2p-Processed\SL1865-Ai203\05142024\SingleP\20PixelFromEdgeExc\';
DataFolder=[ProcessFolder 'Data\'];
MP=load([ProcessFolder 'SLMIncludedIndFromIscell.mat']);
% SingP=[DataFolder 'SingleP\GPL.gpl']
SinglePSTHFolder=[DataFolder 'SinglePSTH\']
ResultFolder=SinglePSTHFolder;
mkdir(SinglePSTHFolder);

% SingPZ=[0 0 50 50 50 100 100]
% BinFile=dir([DataFolder '*TSeries*GPoint*.bin'])
% % SingPZ=[0 0 0 0 50 100 100 100 100 100]
% % SingPZ=[0]
% BinTable = PointLaser_files(BinFile);
TiffTable = ExpInfoMultiTiffFolder(DataFolder);

SessInfo=readtable([ProcessFolder 'SessionInfo.xlsx']);
SessInfo=removevars(SessInfo,{'TotalRepetition'});
SessInfo=removevars(SessInfo,{'Power'});

SessInfo.Properties.VariableNames{'Session'} = 'FileID';
% SessInfo.Properties.VariableNames{'Power'} = 'Laser';
% SessInfo.Properties.VariableNames{'TotalRepetition'} = 'totalRepetitions';
% SessInfo.Properties.VariableNames{'Session'} = 'FileID';


DataList = outerjoin(TiffTable, SessInfo, 'Keys', 'FileID', 'MergeKeys', true);

% DataList = outerjoin(TiffTable, BinTable, 'Keys', 'FileID', 'MergeKeys', true);
[NeuronPos3D,NeuronPos3DRaw,CaData,CaDataPlane,stat,yaml]=Suite2pMultiPlaneROIToXYZ(DataFolder,DataList.FileID(1));
[cellIDMap,CellPixCount,MedCenter,cellBoundary]=Suite2pCellIDMapFromStat(CaData.statCell,[yaml.FOVsize_PX yaml.FOVsize_PY]);
% [cellIDMap,CellPixCount,MedCenter,cellBoundary]=Suite2pCellIDMapFromStat(stat,FovSize)
% cellBoundary = CellIDMap2Boundary(cellIDMap);
yaml.Zdepth_ETL=unique(round(yaml.Zdepth_ETL));
%% In this recording file, I forgot take single image to identify the GPL Position. 
% gplFile='F:\LuSLMOnlineTest\SL0242-Ai203\08202024\Allgpl.gpl';
resultPaths = findAllFiles([ProcessFolder], 'All.gpl');

% resultPaths = findAllFiles([ProcessFolder 'AllIncluded\'], 'GPL.gpl');
% gplFile=resultPaths{1};
% tbl=gpl2Table(gplFile);
% % XYPos=gplXYtoPixel(tbl,yaml);
% [MPPos,Z]=gplXYtoPixel(tbl,yaml);
% MPName=tbl.Name;
% MP3D=[MPPos yaml.Zdepth_ETL(Z(:,2))'];

if exist([ProcessFolder 'SLMIncludedIndFromIscell.mat'])
    MPtemp=load([ProcessFolder 'SLMIncludedIndFromIscell.mat']);
    MP3D=MPtemp.Pos3DNeed;
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
suite2pPath = findAllFolders(DataFolder, 'suite2p');
if length(suite2pPath)==1
   CombinedSuite2p= findAllFiles([suite2pPath{1} 'combined\'], 'Fall.mat');
end
if length(CombinedSuite2p)==1
Fall=load(CombinedSuite2p{1});
cellInfo=Suite2pCellInfo(Fall);
end


if sum(DataList.totalRepetitions)~=size(Fall.F,2)
   disp('Warning!!! Size of Frames from Suite2p does NOT match tiff nums in Tiff Folders');
end
% SinglePxyz=[];
% SinglePxyz=Pos3DNeed;
Laser=unique(DataList.Laser);
Point=unique(DataList.Point);
Laser(isnan(Laser))=[];
Point(isnan(Point))=[];
% % SinglePxyz=Pos3DNeed(Point,:);
SLMframe=median(DataList.PostSLMframe)-1;
% DataList.PostSLMframe(1);

TrialRepTh=4;
PreImgN=10;
% PostImgN=min(90,median(DataList.totalRepetitions));
PostImgN=median(DataList.totalRepetitions);

PostSLMFrameN=3;
PreSLMFrameN=20;


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
NonSLMInd=find(SessInfoNeed.Laser==0&SessInfoNeed.Point>0);
SLMInd=find(SessInfoNeed.Laser(:,1)>1&SessInfoNeed.Point>0);     %% Point  = 0 refers no SLM, < 0 refers to Group Stimuli. 
MarkPoint=unique(SessInfoNeed.Point(SLMInd));

SLM3D=MP3D(Point,:);
SLMName=MPName(Point);



colorCell = jet(length(cellBoundary));
colorCell = colorCell(randperm(length(cellBoundary)),:);


MeanImgClim=[0 1];


figure;
% Img=AmpNormalize(permute(double(CaData.PlaneMeanImg),[2 1 3]),[1 99]);
Img=AmpNormalizeDim(permute(double(CaData.PlaneMeanImg),[2 1 3]),3,[1 99]);
MultiMatrix3DPlotZ(Img,yaml.Zdepth_ETL,0.9);
caxis(MeanImgClim);
Radius=4;
colormap(gray);
set(gca,'xlim',[0 512],'ylim',[0 512],'zlim',yaml.Zdepth_ETL([1 end]),'View',[64 24],'zDir','reverse');
 % plotCellCenter3D(SinglePxyz(:,:,iPoint), Radius, [0 1 0],1.5);
% plotCellCenter3D(NeuronPos3D, Radius, [0 1 0],1);
plotCellCenter3D(SLM3D, Radius, [0 1 0],1);
labelCellCenter(SLM3D, SLMName,[0 1 0])
plotCellBoundary3D(cellBoundary, NeuronPos3D(:,3),colorCell,0.5)
% labelCellCenter(NeuronPos3D, 1:size(NeuronPos3D,1),colorCell)
papersizePX=[0 0 16 length(yaml.Zdepth_ETL)*12];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
saveas(gcf,[ResultFolder '3DCell'],'png');
saveas(gcf,[ResultFolder '3DCell'],'fig');


CellID=1:size(NeuronPos3D,1);

figure;
% labelCellCenter(NeuronPos3D(I,[2 1]), CellID(I),colorCell(I,:));
MultiPlanes2DShow(Img, [], NeuronPos3D, CellID,yaml.Zdepth_ETL, colorCell, MeanImgClim)
papersizePX=[0 0 length(yaml.Zdepth_ETL)*11 10];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
saveas(gcf,[ResultFolder '2DCell'],'png');
saveas(gcf,[ResultFolder '2DCell'],'fig');


%% 
DistTh=10;
[SLMtarget,SLMtargetCellDist]=SLMtargetMatchCell(SLM3D,NeuronPos3D,DistTh)


LaserG=Laser;
% RepG=unique(SessInfoNeed.Repetitions)
% PmtG=unique(SessInfoNeed.PMTLevel)
deltaFoF=F2deltaFoF(Fall.F,Fall.Fneu,Fall.ops.fs);


NData={deltaFoF'};
Nlabel={'DeltaF'};

planeG=unique(cellInfo.iplane);
for iplane=1:length(planeG)
    PlaneC(iplane)=max(find(cellInfo.iplane==planeG(iplane)));
end

iscell=find(Fall.iscell(:,1)==1);
ClimScale=[-5 5;-40 40]
ClimScale=[-1 1;-1 1]
ClimScale=[-2 2;-5 5]
ClimScale=[-0.3 0.3;-0.3 0.3]

ClimScale=[-5 5;-15 15]
   P.xLeft=0.1;        %%%%%%Left Margin
   P.xRight=0.1;       %%%%%%Right Margin
   P.yTop=0.02;         %%%%%%Top Margin
   P.yBottom=0.1;      %%%%%%Bottom Margin
   P.xInt=0.01;         %%%%%%Width-interval between subplots
   P.yInt=0.03;         %%%%%%Height-interval between subplots

PlaneNumStart=[1 PlaneC(1:end-1)+1]
jet=colormap("jet");
colorLaser=jet(1:size(jet,1)/length(LaserG):size(jet,1),:);
close all

%% Ave SLM targets illustration
close all
ResultFolderCell=[ResultFolder 'SLMtarget\'];
mkdir(ResultFolderCell);
% for iCell=1:size(cellInfo,1)

P.xLeft=0.1;        %%%%%%Left Margin
P.xRight=0.1;       %%%%%%Right Margin
P.yTop=0.02;         %%%%%%Top Margin
P.yBottom=0.1;      %%%%%%Bottom Margin
P.xInt=0.01;         %%%%%%Width-interval between subplots
P.yInt=0.03;         %%%%%%Height-interval between subplots
% I1=intersect(NonSLMInd,SingleRep);
close all

PostSLMFrameN=3;
PreSLMFrameN=30;

SLMTable = [];
% SLMframe=SLMframe(1);
LaserG=setdiff(LaserG,0)
for jCell=1:length(SLMtarget)

    iCell=SLMtarget(jCell);
    if iCell<0
       SLMTable(end+1,:) =  [Point(jCell) 0 0 0 1 0 PreSLMFrameN PostSLMFrameN 0 0 0];
       continue;
    end
    figure;
% for iRep=2:2
    for iData=1:length(NData)
        % TempData=zscore(double(NData{iData}(iscell,:)),0,2);
        TempData=double(NData{iData}(iscell(iCell),:));
        TempData=AmpNormalize(TempData,[0 100]);

            for iLaser=1:length(LaserG)
            % for iLaser=1:8
                if LaserG(iLaser)==0
                   I2=find(SessInfoNeed.Laser(:,1)==LaserG(iLaser));      
                else
                   I2=find(SessInfoNeed.Laser(:,1)==LaserG(iLaser)&SessInfoNeed.Point==MarkPoint(jCell));

                end

                I3=I2;
                if length(I3)>=TrialRepTh
                   preSLM=[];
                   postSLM=[];
                   for iSess = 1:length(I3)
                       for iS=1:length(SLMframe)
                       TempI1=SessInfoNeed.FirstFrame(I3(iSess))+SLMframe(iS)-1;
                       PostStim=TempI1:TempI1+PostSLMFrameN-1;
                       PreStim=TempI1-PreSLMFrameN:TempI1-1;
                       % preSLM(iSess)=mean(TempData(PreStim));
                       % postSLM(iSess)=mean(TempData(PostStim));
                       % SLMRatio=mean((postSLM-preSLM)./preSLM);

                       preSLM=[preSLM;TempData(PreStim)'];
                       postSLM=[postSLM;TempData(PostStim)'];
                       end
                   end
                   % SLMRatio=(postSLM-preSLM)./preSLM;
                   % SLMRatio=mean((postSLM-preSLM)./preSLM);
                   % [~,p,~,t]=ttest(SLMRatio,0,'Tail','right')
                   [~,p,~,t]=ttest2(postSLM,preSLM,'Tail','right');
                   ChangePerc=100*(mean(postSLM)-mean(preSLM))/mean(preSLM); 
                   SLMTable(end+1,:) =  [Point(jCell) iCell LaserG(iLaser) ChangePerc p length(I3) PreSLMFrameN PostSLMFrameN t.tstat t.df t.sd ];
                end


                if ~isempty(I3)
                   tempPSTH=[];

                   for iSess = 1:length(I3)
                       TempI=SessInfoNeed.PreFrame(I3(iSess)):SessInfoNeed.PostFrame(I3(iSess));
                       tempPSTH(:,iSess)=TempData(TempI);
                   end
                   PSTHLaser(:,iLaser)=squeeze(mean(tempPSTH,2));
                   % error_area(1:size(tempPSTH,1),mean(tempPSTH,2),std(tempPSTH,0,2),colorLaser(iLaser,:),0.5)
                   subplotLU(1,length(NData),1,iData,P);hold on
                   BaseLine=repmat(mean(tempPSTH(1:PreImgN,:),1),size(tempPSTH,1),1);
                   tempPSTH=tempPSTH-BaseLine;


                   error_area(1:size(tempPSTH,1),mean(tempPSTH,2),ste(tempPSTH')',colorLaser(iLaser,:),0.4,'-',0.5);

                   set(gca,'ylim',[-0.02 0.12]);
                   hold on;
                   for jf=1:length(SLMframe)
                   plot(PreImgN+[SLMframe(jf) SLMframe(jf)],[-0.1 0.3],'k:');
                   end


                   % plot(PSTHLaser(:,iLaser),'color',colorLaser(iLaser,:));
                   % set(gca,'xtick',[0 PreImgN+0.5 PreImgN+0.5+PostImgN],'xticklabel',{'Start' num2str(PreImgN) num2str(PreImgN+PostImgN)})
                   % set(gca,'tickdir','out','clim',ClimScale(iData,:))
                   % if iRep==3
                   set(gca,'xlim',[0 PreImgN+0.5+PostImgN] ,'xtick',[0 PreImgN+0.5 PreImgN+0.5+PostImgN],'xticklabel',{['-' num2str(PreImgN)] '0' num2str(PostImgN)})
                   xlabel(Nlabel{iData})
                   % else
                   % set(gca,'xlim',[0 PreImgN+0.5+PostImgN] ,'xtick',[0 PreImgN+0.5 PreImgN+0.5+PostImgN],'xticklabel',[])
                   % 
                   % end
                   set(gca,'tickdir','out')
                   set(gca,'ylim',[-0.2 0.6],'ytick',[-0.2:0.2:0.6])
                   % if iData==1
                   %    ylabel(['Frame #' num2str(RepG(iRep))]);
                   % end
                   title(['MP' num2str(num2str(MarkPoint(jCell))) 'Cell' num2str(iCell)]);
                   % colormap(ColorPN3)

                end

            end
     colormap(colorLaser)    
     b = colorbar;
     set(b,'position',[0.92 0.2 0.01 0.5]);
     ylabel(b, 'PV power');
     set(b,'ytick',[1:length(LaserG)]/length(LaserG),'yticklabel',LaserG);

    end
     papersizePX=[0 0 length(NData)*8 5 ];
     set(gcf, 'PaperUnits', 'centimeters');
     set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
     saveas(gcf,[ResultFolderCell 'SLMTargetMarkPoint' num2str(MarkPoint(jCell))],'png');
     % close all


end

SLMTable=array2table(SLMTable,'VariableNames',{'MarkPointID','NeuronID','SLMpower','PercChange','p_value','NumTrials','PreSLMframes','PostSLMframes','t','df','sd'})
ZeroPower=find(SLMTable.SLMpower==0);
SLMTable(ZeroPower,:)=[];


writetable(SLMTable,[ResultFolderCell 'SLMtable.xlsx']);
% writetable(SLMTable,[ResultFolderCell 'SLMtable.csv']);
SLMresponse=SLMTable(SLMTable.PercChange>0&SLMTable.p_value<0.05,:);

clear WriteStr
WriteStr{1}=[num2str(length(unique(SLMresponse.MarkPointID))) ' out of ' num2str(length(unique(SLMTable.MarkPointID))) ' SLM targets respond' ' from ' ProcessFolder]
% xlswrite([ResultFolderCell 'SLMtable.xlsx'],WriteStr,'Conclusion');
GlobalResult=['F:\GlobalSLMtargets\SLMtargets.txt'];
FID=fopen(GlobalResult,'a');
fprintf(FID, '%s\r\n', WriteStr{1});
fclose(FID);
GlobalTable=['F:\GlobalSLMtargets\SLMtable.xlsx'];
updateTable(GlobalTable, SLMTable)



