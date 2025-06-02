
function FOVUpdate=BatchSub_AllDataMotionCorrect_PostTiffRemove(FOV,TSeriesBrukerTBL,MotionTh,RemoveFrame)

   % RemoveFrame=2; %%2 repetitions right together with MarkPoint would be removed.

WorkFolder=FOV.WorkFolder;
ConfigFolder='C:\Users\zhangl33\Projects\Project-SLMonlineControl\config\';
ConfigFile='CurrentSLMsetting.yml';%<----------------------------------------------------------------------------------Edit, configuration file
[~,~,~,CaData,CaDataPlane,stat,yaml,confSet]=ROIToXYZ(WorkFolder,ConfigFile);
ProcessFolder = Get_ExpDataFolder(WorkFolder, 'SpeedStimEdgeExc', {'Data','AllIncluded','DataSum','.gpl','.xml'})

PowerTestTiffNum=FOV.confSet.Ziteration*FOV.confSet.ZRepetition*length(FOV.confSet.ETL);
GroupFunTiffNum=sum(TSeriesBrukerTBL{1}.Reps)*length(FOV.confSet.ETL);


SLMPosInfo=load([ProcessFolder 'SLMFunGroup.mat']);
SLMTestInfo=load([ProcessFolder 'SLMIncludedIndFromIscell.mat']);
confSet=SLMPosInfo.confSetFinal;
Zdepth=confSet.scan_Z+confSet.ETL
DataFolder=[ProcessFolder 'Data\'];


ExFolder=[DataFolder 'ExcludeFolder\'];
mkdir(ExFolder);

% % 
% % 
[FileGenerateInfo,fileList, fileIDs,tiffNum] = getExpInfoFiles_NonMat(DataFolder);
PowerTestfileI=find(tiffNum==PowerTestTiffNum)
GroupFunfileI=find(tiffNum==GroupFunTiffNum)

PowerTestSLMtbl=[];
PowerTestSLMNum=[];
for iFile=1:length(PowerTestfileI)
    % [OutTBLTemp,XYTrials,PowerWeight]=MPSeqFolder_GroupTargets(Folder,SLMTestInfo)
    OutTBLTemp=MPSeqFolder_1TargetXNon([DataFolder fileList{PowerTestfileI(iFile)} '\'],[confSet.SLM_Pixels_Y;confSet.SLM_Pixels_X],SLMTestInfo.Pos3Dneed);
    OutTBLTemp.FileID = fileIDs(PowerTestfileI(iFile)) * ones(size(OutTBLTemp, 1), 1);
    PowerTestSLMtbl=[PowerTestSLMtbl;OutTBLTemp];
    PowerTestSLMNum(iFile,:)=[fileIDs(PowerTestfileI(iFile)) size(OutTBLTemp,1)];
end
% % 
% % 
Pos3DFun=SLMPosInfo.FinalPos3D;
% FunScore=SLMPosInfo.FinalFunScore;
Group=SLMPosInfo.Group;
for iGroup=1:length(Group)
    Pos3DGroup{iGroup}=Pos3DFun(Group(iGroup).Indices,:);
end


FunSLMtbl=[];
FunSLMNum=[];
for iFile=1:length(GroupFunfileI)
    % [OutTBLTemp,XYTrials,PowerWeight]=MPSeqFolder_GroupTargets(Folder,SLMTestInfo)
    OutTBLTemp=MPSeqFolder_GroupTargets([DataFolder fileList{GroupFunfileI(iFile)} '\'],[confSet.SLM_Pixels_Y;confSet.SLM_Pixels_X],Pos3DGroup);
    OutTBLTemp.FileID = fileIDs(GroupFunfileI(iFile)) * ones(size(OutTBLTemp, 1), 1);
    FunSLMtbl=[FunSLMtbl;OutTBLTemp];
    FunSLMNum(iFile)=size(OutTBLTemp,1);
end



% % 
% % 
InvalidI1=PowerTestfileI(PowerTestSLMNum(:,2)<(confSet.Ziteration-1));
InvalidI2=PowerTestfileI(FunSLMNum(:,2)<sum(TSeriesBrukerTBL{1}.SynMP));
InvalidI=union(InvalidI1,InvalidI2);
copyInitialRecordedFolders(fileList(InvalidI), ExFolder, DataFolder);
DelFolders(fileList(InvalidI), DataFolder);
% % 
% % 
FunSLMtbl = MatchOutTBLAll_TSeriesBruker(FunSLMtbl, TSeriesBrukerTBL);

T = outerjoin(PowerTestSLMtbl, FunSLMtbl, 'MergeKeys', true);
[~,I1]=sort(T.FileID);
T=T(I1,:);

FileInfo=table(fileIDs(:),tiffNum(:), fileList(:),'VariableNames',{'FileID', 'tiffNum','FileKey'});

FileInfo = innerjoin(FileInfo, FOV.MatFile, "Keys", "FileID");

TT = outerjoin(T, FileInfo, "Keys", "FileID", "Type", "right", "MergeKeys", true);



OutTBLAll=TT; 
OutTBLAll.AwakeState(TT.TSeriesInd<=5)=1;    %%The 1st half 5 Tseries is designed for awake state.
OutTBLAll.AwakeState(TT.TSeriesInd>=6)=2;    %%The 2nd half 5 Tseries is designed for anesia state.


[FileID,I1]=unique(OutTBLAll.FileID);
InvalidI=OutTBLAll.motionMed(I1)>MotionTh&OutTBLAll.AwakeState(I1)~=2;         %%Remove motion shift large recording except for anesia sessions.
AllFileList=OutTBLAll.FileKey(I1);
copyInitialRecordedFolders(AllFileList(InvalidI), ExFolder, DataFolder);
DelFolders(AllFileList(InvalidI), DataFolder);

OutTBLAll(ismember(OutTBLAll.FileKey,AllFileList(InvalidI)),:)=[];


[FileID,~]=unique(OutTBLAll.FileID);

[FOVUpdate.FileGenerateInfo,FOVUpdate.fileList, FOVUpdate.fileIDs,FOVUpdate.tiffNum] = getExpInfoFiles_NonMat(DataFolder);

% % 
[~,~, fileIDCurrent,tiffNum] = getExpInfoFiles_NonMat(DataFolder);
if isempty(setdiff(fileIDCurrent,FileID))&&isempty(setdiff(FileID,fileIDCurrent))
   disp('FileID in tiff folder and OutTBLAll Table match, continue to delete post SLM tiff files for all folders');
   [TiffTable, RemoveList] = RemoveMPsynTiffFolder(DataFolder,RemoveFrame);

    [FileGenerateInfo,fileList, fileIDs,tiffNum] = getExpInfoFiles_NonMat(DataFolder);
    % [FOVUpdate.FileGenerateInfo,FOVUpdate.fileList, FOVUpdate.fileIDs,FOVUpdate.tiffNum] = getExpInfoFiles_NonMat(DataFolder);

    PostTiffTable=table(fileIDs(:),tiffNum(:),'VariableNames',{'FileID','Suite2pTiffNum'})
    Suite2pTable=outerjoin(OutTBLAll, PostTiffTable, "Keys", "FileID", "Type", "right", "MergeKeys", true); 
    save([DataFolder 'TableForSuite2p.mat'],'Suite2pTable','SLMPosInfo','SLMTestInfo','RemoveFrame');

    % FOVUpdate=FOV;
    FOVUpdate.DataFolder=DataFolder;
    FOVUpdate.Suite2pTable=Suite2pTable;
    FOVUpdate.FileGenerateInfo=FileGenerateInfo;
    FOVUpdate.fileList=fileList;
    FOVUpdate.fileIDs=fileIDs;
    FOVUpdate.tiffNum=tiffNum;
    FOVUpdate.Suite2pTable=Suite2pTable;

    [~,i1]=unique(Suite2pTable.FileID);
    subT1=Suite2pTable(i1,:);

    FOVUpdate.subT1=subT1;

    disp(['Total of ' num2str(sum(subT1.Suite2pTiffNum)) ' tif files required to processed in suite2p']);
else
    disp('FileID in tiff folder and OutTBLAll Table do NOT match, check!');
end
% % 
