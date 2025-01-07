%% Initialize, run this part before mannual correction of Suite2p processed data for saving time
function [rSpeed,rStim,nFrame,fSpeedOri]=FunScore_Cal(LoadPath,fs)
% LoadPath='\\nimhlabstore1.nimh.nih.gov\UFNC\FNC3Z\RawRecording\SL0777-Ai203\10292024\';

% ConfigFolder='C:\Users\User\Project-SLMonlineControl\config\';
ConfigFolder='C:\Users\zhangl33\Projects\Project-SLMonlineControl\config\';

ConfigFile='SLMsetting.yml';
confSet = ReadYaml([ConfigFolder '\' ConfigFile]);
confSet.save_path0=LoadPath;
fs=6.9;



% This part cost 3 minutes.
if ~exist('fSpeed')
    [fSpeed,fStim,timeStampCa_Plane]=PV_VolExtract(confSet);
    fSpeedOri=fSpeed;
    fSpeed=AmpNormalizeRow(double(fSpeed)',[0 100])';
    fStim=AmpNormalizeRow(double(fStim)',[0 100])';
end

[Pos3D,Pos3DRaw,CaData,CaDataPlane,stat,yaml]=ROIToXYZ_FromPath(LoadPath);
PlaneZ=yaml.Zdepth_ETL+yaml.scan_Z(1);


IndexNeed=1:1:size(Pos3D,1);
% XYZtoMarkPoint(SavePathAllPoint,Pos3D,IndexNeed,yaml,confSet,CaData.statCell);
numPlanes=length(confSet.ETL);
iscell=find(CaData.iscell(:,1)==1);

% Exlude cells near the edge of the FOV as SLM targets
numPoint=size(Pos3D,1);
% XYrange=[SLMrangePix;yaml.SLM_Pixels_Y-SLMrangePix];  %%Cell locates close to edge of the view, were not considered as SML targets.
% OutRange=find(Pos3D(:,1)<XYrange(1)|Pos3D(:,2)<XYrange(1)|Pos3D(:,1)>XYrange(2)|Pos3D(:,2)>XYrange(2));
% CenterCloseI=setdiff(1:numPoint,OutRange);

[rSpeed,rStim]=FunFilterCell(CaData,iscell,fSpeed,fStim,fs);
nFrame=size(CaData.F,2);


%% After Suite2p processeing is done
% Intially all cells were dectected by suite2p were considered as SLM targets
% SLMrangePix=20; %Pixel number close to FOV is excluded
% Exlude cells near the edge of the FOV as SLM targets
% step3_SubStep1_FOVedgeExcludedROIMP;    %% CenterCloseI is generated at this step


save([LoadPath 'FunScoreBeh.mat'],'fStim','fSpeed','rSpeed','rStim');


%% step3_generateMP from steps above without Functional Filter Exlude cells near the edge of the FOV as SLM targets
% SLMrangePix=50; %Pixel number close to FOV is excluded
% numPoint=size(Pos3D,1);
% XYrange=[SLMrangePix;yaml.SLM_Pixels_Y-SLMrangePix]  %%Cell locates close to edge of the view, were not considered as SML targets.
% OutRange=find(Pos3D(:,1)<XYrange(1)|Pos3D(:,2)<XYrange(1)|Pos3D(:,1)>XYrange(2)|Pos3D(:,2)>XYrange(2));
% CenterCloseI=setdiff(1:numPoint,OutRange);
% SavePathExc=[SavePath num2str(SLMrangePix) 'PixelFromEdgeExc\']
% mkdir(SavePathExc)
% XYZtoMarkPoint_NT_PairGplXml(SavePathExc, Pos3D, CenterCloseI, NonTargets, IndexNonTargetTrial, yaml, confSet,CaData.statCell);



%% Including top cells highly correlated associated with speed.
% TopCellN=13;  
% SavePathStimSpeed=[SavePath 'Top' num2str(TopCellN) 'SpeedStimEdgeExc\'];
% mkdir(SavePathStimSpeed);
% 
% step3_SubStep4_FuncFilterMPs;




%% After Mannual Adjust Generated Functional Filtered MarkPoints in PV
% PostGPLadjusted(SavePathStimSpeed,IncludeCellFunFilter, FunScore, yaml, confSet,CaData.statCell,NonTargetPath,NonTargets, IndexNonTargetTrial);
%%




% figure;
% plot(rSpeed(CenterCloseI),'r-');
% hold on;
% plot(rStim(CenterCloseI),'g-');
% plot(rCenterISpeed(1:TopCellN),0.3,'ro')
% plot(rCenterIStim(1:TopCellN),0.3,'g*')
% 
% 
% 
% figure;
% for iTop=1:TopCellN
% % subplot(TopCellN,1,i)
% tempData=NeuroData(:,iscell(CenterCloseI(rCenterIStim(iTop))));
% tempData=SmoothDec(tempData,1);
% plot(tempData+iTop);
% hold on;
% plot(mean(fStim,2)*(TopCellN+1),'r');
% end
% 
% 
% 
% figure;
% for iTop=1:TopCellN
% % subplot(TopCellN,1,i)
% tempData=NeuroData(:,iscell(CenterCloseI(rCenterISpeed(iTop))));
% tempData=SmoothDec(tempData,1);
% plot(tempData+iTop);
% hold on;
% plot(mean(fSpeed,2)*(TopCellN+1),'r');
% end
% 
% 
% subplot(2,1,1)
% imagesc(NeuroData(:,iscell(CenterCloseI))');
% subplot(2,1,2)
% plot(mean(fSpeed,2));hold on;
% 
% 
% figure;
% for i=1:length(CenterCloseI)
%     subplot(length(CenterCloseI),1,i)
% plot(NeuroData(:,iscell(CenterCloseI(i))));hold on;
% plot(mean(fStim,2)+5,'r-');
% end


% figure;
% plot(rSpeed,'r');
% hold on;
% plot(rStim,'g');
% hold on;
% plot(rFspks,'b')

% figure;
% plot(-log10(pFspks),'b')
% [~,r1]=sort(rFspks(CenterCloseI),'descend')  %%Cell locates close to edge of the view, were not considered as SML targets.
% TopCellN=8;         %%   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<Top X cells highly associated with speed 
% TopSpeedCellI=CenterCloseI(r1(1:TopCellN));
% SavePathSpeed=[SavePath 'Top' num2str(TopCellN) 'SpeedEdgeExc\']
% mkdir(SavePathSpeed)
% XYZtoMarkPoint_NT(SavePathSpeed,Pos3D,TopSpeedCellI,NonTargets, IndexNonTargetTrial, yaml,confSet)
% 





