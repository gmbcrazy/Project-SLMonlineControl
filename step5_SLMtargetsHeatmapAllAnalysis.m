

clear all


load('C:\Users\zhangl33\Projects\Project-SLMonlineControl\subfun\Color\colorMapPN3.mat')
% Param.TrialRepTh=4;
Param.PreImgN=0;
% PostImgN=min(90,median(DataList.totalRepetitions));

Param.PreSLMFrameN=30; %Frame num. before SLM to calculate baseline level
Param.PostSLMFrameN=3; %Frame num. after SLM to calculate response level
Param.TrialRepTh=4; %%Minimal num. of trials to quantify SLM response
Param.DistTh=10;%%Maximal distance by pixels considered between target cell and SLM target center.
Param.PixFromMedCenter=20;
Param.Clim=[-300 300];
Param.ColorMap=colorMapPN1;



TargetZ=[];
SessinInfo=[];
GlobalSavePath='\\nimhlabstore1.nimh.nih.gov\UFNC\FNC2\Zhang\SLM\SummarySLMResponse\AllCurentSLMtarget';


ProcessFolder='F:\LuSLMOnlineTest\SL0543-Ai203\09102024\SingleP\Top8SpeedStimEdgeExc\';
SessinInfo=[];
SLMTable=SLMtargetAnalysis(ProcessFolder, Param,GlobalSavePath,SessinInfo, TargetZ);


ProcessFolder='F:\LuSLMOnlineTest\SL0242-Ai203\09042024\SingleP\Top5SpeedStimEdgeExc\';
SessinInfo=[];
[HMap,RoiMap]=SLMtargetHeatMapAnalysis(ProcessFolder, Param,GlobalSavePath,SessinInfo, TargetZ);

ProcessFolder='F:\LuSLMOnlineTest\SL0242-Ai203\08292024\SingleP\Top5SpeedStimEdgeExc\';
[HMap,RoiMap]=SLMtargetHeatMapAnalysis(ProcessFolder, Param,GlobalSavePath,SessinInfo, TargetZ);


threshold_percentage=0.2;
thNum=15;
max_distance=12;
min_region_size=10;
PeakTh=200;
min_merged_region_size=20;
figure;
for iCell=1:length(RoiMap)
    if isempty(RoiMap{iCell})
       continue;
    end
    resultTemp=[];
    for jLaser=1:length(RoiMap{iCell})
        subplotLU(2,length(RoiMap),jLaser,iCell)
        tempData=RoiMap{iCell}{jLaser};
        % Cell=FieldFind2D_Lu(tempData,RateTH,thNum);
        % [result, regions] = check_high_value_center(tempData, threshold_percentage, max_distance, min_region_size, PeakTh)
        % [result, mergedRegion, contourPixels] = check_high_value_center(tempData, threshold_percentage, max_distance, min_region_size, PeakTh);
        [resultTemp(jLaser), mergedRegion, contourPixels] = check_high_value_center(tempData, threshold_percentage, max_distance, min_region_size, PeakTh, min_merged_region_size);
        imagesc(tempData);hold on;
        set(gca,'clim',[-300 300])
        colormap(colorMapPN1)
        % for iField=1:length(regions)
        if ~isempty(mergedRegion)
            plot(contourPixels(:,2),contourPixels(:,1),'g.');
        end
        % end

    end
    result(iCell)=max(resultTemp);
end




Cell=FieldFind2D_Lu(mapO,RateTH,thNum)



ProcessFolder='D:\LuSLMOnlineTest\SL0340\05302024\Ch2\TSeries\PowerLevel1\';
SessinInfo=1;
SLMTable=SLMtargetAnalysis(ProcessFolder, Param,GlobalSavePath,SessinInfo, TargetZ);


ProcessFolder='\\nimhlabstore1.nimh.nih.gov\UFNC\FNC3Z\Sutie2p-Processed\SL1865-Ai203\05142024\SingleP\20PixelFromEdgeExc\';
SessinInfo=1;
SLMTable=SLMtargetAnalysis(ProcessFolder, Param,GlobalSavePath,SessinInfo, TargetZ);



ProcessFolder='F:\LuSLMOnlineTest\SL0541-Emx1Ai96\09122024\SingleP\Top8SpeedStimEdgeExc\';
SessinInfo=[];
SLMTable=SLMtargetAnalysis(ProcessFolder, Param,GlobalSavePath,SessinInfo, TargetZ);

ProcessFolder='F:\LuSLMOnlineTest\SL0541-Emx1Ai96\09172024\SingleP\Top10SpeedStimEdgeExcTest\';
SessinInfo=[];
SLMTable=SLMtargetAnalysis(ProcessFolder, Param,GlobalSavePath,SessinInfo, TargetZ);



ProcessFolder='F:\LuSLMOnlineTest\SL0702-G7-Ai203\09192024\SingleP\Top15SpeedStimEdgeExc\';
SessinInfo=[];
SLMTable=SLMtargetAnalysis(ProcessFolder, Param,GlobalSavePath,SessinInfo, TargetZ);
