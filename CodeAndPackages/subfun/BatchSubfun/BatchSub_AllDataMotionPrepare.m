function [MatFile,FileGenerateInfo,fileList, fileIDs,tiffNum,confSet]=BatchSub_AllDataMotionPrepare(WorkFolder,TSeriesBrukerTBL)

% Copy 1st two Spontanous Behavior Sessions to later PowerTest and SLMgroup Data folder
% Delete Recording sessions with wrong Tiff num (Does not match spontnous recording, PowerTest, SLMgroup experiment)
% Also return the motion pixel shifts from all files.

% WorkFolder='E:\LuSLMOnlineTest\SL0864\04162025\';
% ConfigFolder='C:\Users\zhangl33\Projects\Project-SLMonlineControl\config\';
ConfigFile='CurrentSLMsetting.yml';%<----------------------------------------------------------------------------------Edit, configuration file
[~,~,~,CaData,CaDataPlane,stat,yaml,confSet]=ROIToXYZ(WorkFolder,ConfigFile);
ProcessFolder = Get_ExpDataFolder(WorkFolder, 'SpeedStimEdgeExc', {'Data','AllIncluded','DataSum','.gpl','.xml'});



SLMPosInfo=load([ProcessFolder 'SLMFunGroup.mat']);
SLMTestInfo=load([ProcessFolder 'SLMIncludedIndFromIscell.mat']);
confSet=SLMPosInfo.confSetFinal;
Zdepth=confSet.scan_Z+confSet.ETL
DataFolder=[ProcessFolder 'Data\'];

%%Get Initial spontanous recording motion shifts, and move tiff files folder to final data folder for further processing
[FileGenerateInfo,InitialfileList, InitialfileID] = getExpInfoFiles_NonMat(WorkFolder);
[InitialPixShiftFile, Files] = PixShiftLoad(WorkFolder);
copyInitialRecordedFolders(InitialfileList, DataFolder,WorkFolder);








%%Exclude sessions with non-correct num. of tiff
[FileGenerateInfo,fileList, fileIDs,tiffNum] = getExpInfoFiles_NonMat(DataFolder);
PowerTestTiffNum=confSet.Ziteration*confSet.ZRepetition*length(confSet.ETL);
GroupFunTiffNum=sum(TSeriesBrukerTBL{1}.Reps)*length(confSet.ETL);


ValidTiffNum=unique(tiffNum(ismember(fileIDs,InitialfileID)|tiffNum==PowerTestTiffNum|tiffNum==GroupFunTiffNum))
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
    for i=1:length(I1)
        [PixShiftFile(i), ~] = PixShiftLoad([DataFolder fileList{I1(i)}]);
    end
    MissingFileID=round(fileIDs(I1));
    NewData = table(MissingFileID(:), PixShiftFile(:), ...
    'VariableNames', {'FileID', 'motionMed'});

    MatFile = [NewData;MatFile]; 
end



