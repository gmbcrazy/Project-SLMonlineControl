clear all

BatchSavePath='D:\Project1-LocalProcessing\Step1\';

Suite2pSaveFolderAll='\\nimhlabstore1.nimh.nih.gov\UFNC\FNC3Z\Sutie2p-Processed\GCamP6S-CamKII\';
Suite2pSaveFolderAllLocal='C:\GCamP6S-CamKII\';

rootData='E:\LuSLMOnlineTest\';
AnimalPath=dir([rootData 'SL????']);
ConfigFolder='C:\Users\zhangl33\Projects\Project-SLMonlineControl\config\';

if exist([BatchSavePath 'FOV.mat'])
   load([BatchSavePath 'FOV.mat']);
end



FOVsession={};
for i = 1:length(AnimalPath)
    parentPath = fullfile(rootData, AnimalPath(i).name);
    allItems = dir(parentPath);

    % Filter: folders with 8-digit names only
    isValid = [allItems.isdir] & ...
              ~ismember({allItems.name}, {'.', '..'}) & ...
              cellfun(@(x) ~isempty(regexp(x, '^\d{8}$', 'once')), {allItems.name});

    dateFolders = allItems(isValid);

    % Example: show full paths
    fullPaths = fullfile(parentPath, {dateFolders.name});
    disp(fullPaths);
    FOVsession=[FOVsession;fullPaths(:)];
end


load([ConfigFolder 'PreGenerateTseriesMultiZ\SpontBeh5T_Z11Frame550.mat'],'TSeriesBrukerTBL');
TSeriesBrukerTBL1=TSeriesBrukerTBL;
load([ConfigFolder 'PreGenerateTseriesMultiZ\Anesthesia5T_Z11Frame550.mat'],'TSeriesBrukerTBL');
TSeriesBrukerTBL2=TSeriesBrukerTBL;
clear TSeriesBrukerTBL
TSeriesBrukerTBL=[TSeriesBrukerTBL1 TSeriesBrukerTBL2];
if exist('FOVUpdate')
   FOVExist=FOVUpdate;
   for i=1:length(FOVExist)
       FOVExistName{i}=[rootData animalIDs{i} '\' sessionDates{i}];
   end
else 
   clear FOV;
end

FOV=struct([]);
UpdateIn=[];
clear FOV;
for i=1:length(FOVsession)
    WorkFolder=[FOVsession{i} '\'];
    if exist('FOVExistName')
       [~,~,I1]=intersect(FOVsession{i},FOVExistName)
       if ~isempty(I1)
          % FOV(i)=FOVUpdate(I1);
          UpdateIn=[UpdateIn;i];
       else
          FOV(i).WorkFolder=WorkFolder;
         [FOV(i).MatFile,FOV(i).FileGenerateInfo,FOV(i).fileList, FOV(i).fileIDs,FOV(i).tiffNum,FOV(i).confSet]=BatchSub_AllDataMotionPrepare(WorkFolder,TSeriesBrukerTBL);
  
          % disp(['Check existing processed tiff folder data of ' FOVsession{i} ' and ' FOVExistName]);
       end
       continue;
       % [FOV(i).MatFile,FOV(i).FileGenerateInfo,FOV(i).fileList, FOV(i).fileIDs,FOV(i).tiffNum,FOV(i).confSet]=BatchSub_AllDataMotionPrepare(WorkFolder,TSeriesBrukerTBL);
    else
        [FOV(i).MatFile,FOV(i).FileGenerateInfo,FOV(i).fileList, FOV(i).fileIDs,FOV(i).tiffNum,FOV(i).confSet]=BatchSub_AllDataMotionPrepare(WorkFolder,TSeriesBrukerTBL);

    end
end



PowerTestTiffNum=FOV(i).confSet.Ziteration*FOV(i).confSet.ZRepetition*length(FOV(i).confSet.ETL)
GroupFunTiffNum=sum(TSeriesBrukerTBL{1}.Reps)*length(FOV(i).confSet.ETL)




UniversalMotionTh=6;
figure;
for i=1:length(FOVsession)
    if ~isempty(FOV(i).WorkFolder)
    subplotLU(length(FOVsession),2,i,1);
    % hist(FOV(i).MatFile.motionMed,20);
    histPlotLU(FOV(i).MatFile.motionMed,0:0.5:20,[0.5 0.5 0.5],0.6);
 
    text(8,10,FOVsession{i},'FontSize',8)

    hold on;
    plot([UniversalMotionTh UniversalMotionTh],[0 20],'r:');


    BadInd=find(FOV(i).MatFile.motionMed>UniversalMotionTh);
    subplotLU(length(FOVsession),2,i,2);
    plot(FOV(i).MatFile.FileID,FOV(i).MatFile.motionMed,'o-');
    hold on;
    plot(FOV(i).MatFile.FileID(BadInd), FOV(i).MatFile.motionMed(BadInd),'ro','MarkerFaceColor','r');

    plot([0 length(FOV(i).MatFile.motionMed)],[UniversalMotionTh UniversalMotionTh]);
    end
end
    xlabel('FileID')

 papersizePX=[0 0 20 length(FOVsession)*8];
 set(gcf, 'PaperUnits', 'centimeters');
 set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
 print(gcf, [BatchSavePath date 'MotionCorrectionControl.svg'], '-dsvg', '-painters');
 print(gcf, [BatchSavePath date 'MotionCorrectionControl.tif'], '-dtiffn', '-painters');
 close all




RemoveFrame=2;
iTemp=1;
for i=1:length(FOV)
    if ~isempty(FOV(i).MatFile)
       FOVUpdate(i)=BatchSub_AllDataMotionCorrect_PostTiffRemove(FOV(i),TSeriesBrukerTBL,UniversalMotionTh,RemoveFrame);
    else
       disp([FOVsession(i) ' previously processed']);
       FOVupdate(i)=FOVExist(UpdateIn(iTemp));
       iTemp=iTemp+1;
    end
end


animalpaths={FOVUpdate.DataFolder};

for i = 1:length(animalpaths)
    % Extract animal ID (e.g., SL0886)
    animalMatch = regexp(animalpaths{i}, 'SL\d{4}', 'match');
    if ~isempty(animalMatch)
        animalIDs{i} = animalMatch{1};
    else
        animalIDs{i} = '';
    end

    % Extract date (8-digit number)
    dateMatch = regexp(animalpaths{i}, '\d{8}', 'match');
    if ~isempty(dateMatch)
        sessionDates{i} = dateMatch{1};
    else
        sessionDates{i} = '';
    end

    suite2pFOVPath{i}=[Suite2pSaveFolderAll '\' animalIDs{i} sessionDates{i} '\'];
    mkdir(suite2pFOVPath{i})
    suite2pFOVPathLocal{i}=[Suite2pSaveFolderAllLocal '\' animalIDs{i} sessionDates{i} '\'];
    mkdir(suite2pFOVPathLocal{i})

end





for i=1:length(FOVUpdate)
    writetable(FOVUpdate(i).subT1,[BatchSavePath 'FOV' num2str(i) animalIDs{i} sessionDates{i} '.csv']);
    FOVUpdate(i).CSVName=[BatchSavePath 'FOV' num2str(i) animalIDs{i} sessionDates{i} '.csv'];
end

save([BatchSavePath date 'FOV.mat'],'FOVUpdate','suite2pFOVPath','suite2pFOVPathLocal');
