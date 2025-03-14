clear all

WorkFolder='E:\LuSLMOnlineTest\SL0855-Emx1G6CII-AAV9CAMKII\03062025\';
ConfigFolder='C:\Users\zhangl33\Projects\Project-SLMonlineControl\config\';
ConfigFile='CurrentSLMsetting.yml';%<----------------------------------------------------------------------------------Edit, configuration file
[~,~,~,CaData,CaDataPlane,stat,yaml,confSet]=ROIToXYZ(WorkFolder,ConfigFile);
ProcessFolder = Get_ExpDataFolder(WorkFolder, 'SpeedStimEdgeExc', {'Data','AllIncluded','DataSum','.gpl','.xml'})

motionTh=5;

SLMPosInfo=load([ProcessFolder 'SLMFunGroup.mat']);
SLMTestInfo=load([ProcessFolder 'SLMIncludedIndFromIscell.mat']);
confSet=SLMPosInfo.confSetFinal;
Zdepth=confSet.scan_Z+confSet.ETL
DataFolder=[ProcessFolder 'Data\'];

%%Get Initial spontanous recording motion shifts, and move tiff files folder to final data folder for further processing
[FileGenerateInfo,InitialfileList, InitialfileID] = getExpInfoFiles_NonMat(WorkFolder)
[InitialPixShiftFile, Files] = PixShiftLoad(WorkFolder);
copyInitialRecordedFolders(InitialfileList, DataFolder,WorkFolder);

load([ConfigFolder '\PreGenerateTseriesMultiZ\SpontBeh5T_Z11Frame550.mat','TSeriesBrukerTBL']);
TSeriesBrukerTBL1=TSeriesBrukerTBL;
load([ConfigFolder 'PreGenerateTseriesMultiZ\Anesthesia5T_Z11Frame550.mat','TSeriesBrukerTBL']);
TSeriesBrukerTBL2=TSeriesBrukerTBL;
clear TSeriesBrukerTBL
TSeriesBrukerTBL=[TSeriesBrukerTBL1 TSeriesBrukerTBL2];




%%Exclude sessions with non-correct num. of tiff
[FileGenerateInfo,fileList, fileIDs,tiffNum] = getExpInfoFiles_NonMat(DataFolder);

PowerTestTiffNum=confSet.Ziteration*confSet.ZRepetition*length(confSet.ETL);
GroupFunTiffNum=sum(TSeriesBrukerTBL{1}.Reps)*length(confSet.ETL);


ValidTiffNum=unique(tiffNum(ismember(fileIDs,InitialfileID)|tiffNum==confSet.Ziteration*confSet.ZRepetition*length(confSet.ETL)|tiffNum==sum(TSeriesBrukerTBL{1}.Reps)*length(confSet.ETL)))
% [FileGenerateInfo,fileList, fileIDs,tiffNum] = getExpInfoFiles_NonMat(DataFolder);
ismember(tiffNum,ValidTiffNum);
Invalidfile=~ismember(tiffNum,ValidTiffNum);
ExFolder=[DataFolder 'ExcludeFolder\'];
mkdir(ExFolder);
copyInitialRecordedFolders(fileList(Invalidfile), ExFolder, DataFolder);
DelFolders(fileList(Invalidfile), DataFolder) 

%%Exclude sessions with non-correct num. of tiff
% [FileGenerateInfo,fileList, fileIDs,tiffNum] = getExpInfoFiles_NonMat(DataFolder);
[MatFile, MatExp] = ExtractExp_FromMat(DataFolder);
[FileGenerateInfo,fileList, fileIDs,tiffNum] = getExpInfoFiles_NonMat(DataFolder);


NewData = table(round(InitialfileID(:)), repmat(InitialPixShiftFile, length(InitialfileID), 1), ...
    'VariableNames', {'FileID', 'motionMed'});

MatFile = [NewData;MatFile]; 

%%Check if there is tiff folder where no experimental .mat files is record (due to wrong deletion when recording)
[~,I1]=setdiff(fileIDs,MatFile.FileID);

MissingPixShiftFile=[];
MissingFileID=[];
if ~isempty(I1)
    for i=length(I1)
        [PixShiftFile(i), ~] = PixShiftLoad([DataFolder fileList{I1(i)}]);
    end
    MissingFileID=round(fileIDs(I1));
    NewData = table(MissingFileID(:), PixShiftFile(:), ...
    'VariableNames', {'FileID', 'motionMed'});

    MatFile = [NewData;MatFile]; 
end


Invalidfile=ismember(fileIDs,MatFile.FileID(MatFile.motionMed>motionTh)');
copyInitialRecordedFolders(fileList(Invalidfile), ExFolder, DataFolder);
DelFolders(fileList(Invalidfile), DataFolder);


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


InvalidI1=PowerTestfileI(PowerTestSLMNum(:,2)<(confSet.Ziteration-1));
InvalidI2=PowerTestfileI(FunSLMNum(:,2)<sum(TSeriesBrukerTBL{1}.SynMP));
InvalidI=union(InvalidI1,InvalidI2);
copyInitialRecordedFolders(fileList(InvalidI), ExFolder, DataFolder);
DelFolders(fileList(InvalidI), DataFolder);


FunSLMtbl = MatchOutTBLAll_TSeriesBruker(FunSLMtbl, TSeriesBrukerTBL);

T = outerjoin(PowerTestSLMtbl, FunSLMtbl, 'MergeKeys', true);
[~,I1]=sort(T.FileID);
T=T(I1,:);

FileInfo=table(fileIDs(:),tiffNum(:), fileList(:),'VariableNames',{'FileID', 'tiffNum','FileKey'});

FileInfo = innerjoin(FileInfo, MatFile, "Keys", "FileID");

TT = outerjoin(T, FileInfo, "Keys", "FileID", "Type", "right", "MergeKeys", true);



OutTBLAll=TT; 
OutTBLAll.AwakeState(TT.TSeriesInd<=5)=1;    %%The 1st half 5 Tseries is designed for awake state.
OutTBLAll.AwakeState(TT.TSeriesInd>=6)=2;    %%The 2nd half 5 Tseries is designed for anesia state.


FileID=unique(OutTBLAll.FileID)

[~,~, fileIDCurrent,tiffNum] = getExpInfoFiles_NonMat(DataFolder);
if isempty(setdiff(fileIDCurrent,FileID))&&isempty(setdiff(FileID,fileIDCurrent))
   disp('FileID in tiff folder and OutTBLAll Table match, continue to delete post SLM tiff files for all folders');
   RemoveFrame=2; %%2 repetitions right together with MarkPoint would be removed.
    [TiffTable, RemoveList] = RemoveMPsynTiffFolder(DataFolder,RemoveFrame);

    [FileGenerateInfo,fileList, fileIDs,tiffNum] = getExpInfoFiles_NonMat(DataFolder);
    PostTiffTable=table(fileIDs(:),tiffNum(:),'VariableNames',{'FileID','Suite2pTiffNum'})
    Suite2pTable=outerjoin(OutTBLAll, PostTiffTable, "Keys", "FileID", "Type", "right", "MergeKeys", true); 
    save([DataFolder 'TableForSuite2p.mat'],'Suite2pTable','SLMPosInfo','SLMTestInfo','RemoveFrame');

    [~,i1]=unique(Suite2pTable.FileID);
    subT1=Suite2pTable(i1,:);
    disp(['Total of ' num2str(sum(subT1.Suite2pTiffNum)) ' tif files required to processed in suite2p']);
else
    disp('FileID in tiff folder and OutTBLAll Table do NOT match, check!');
end


