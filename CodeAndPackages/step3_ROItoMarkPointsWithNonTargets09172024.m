%% Initialize, run this part before mannual correction of Suite2p processed data for saving time
clear all
ConfigFolder='C:\Users\User\Project-SLMonlineControl\config\';

ConfigFile='SLMsetting.yml';
confSet = ReadYaml([ConfigFolder '\' ConfigFile]);

% This part cost 3 minutes.
fileID=[1 2];
if ~exist('fSpeed')
    [fSpeed,fStim,timeStampCa_Plane]=PV_VolExtract_MultiFolder(confSet,fileID);
    % [fSpeed,fStim,timeStampCa_Plane]=PV_VolExtract(confSet);
    fSpeed=AmpNormalizeRow(double(fSpeed)',[0 100])';
    fStim=AmpNormalizeRow(double(fStim)',[0 100])';
end


%% After Suite2p processeing is done
% Intially all cells were dectected by suite2p were considered as SLM targets
SLMrangePix=20; %Pixel number close to FOV is excluded
% Exlude cells near the edge of the FOV as SLM targets
step3_SubStep1_FOVedgeExcludedROIMP;    %% CenterCloseI is generated at this step


%% Automatic generate Raw Non-Targets 
step3_SubStep2_RawNonTarget


%% Please do mannual correction to exclude disqualified non-targets in PV, after that, exported all selected targets to SelectedFromRaw.gpl file
step3_SubStep3_AfterMannualSelectedNonTargetsInPV




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
TopCellN=12;  
SavePathStimSpeed=[SavePath 'Top' num2str(TopCellN) 'SpeedStimEdgeExc\'];
mkdir(SavePathStimSpeed);

step3_SubStep4_FuncFilterMPs;




%% After Mannual Adjust Generated Functional Filtered MarkPoints in PV
PostGPLadjusted(SavePathStimSpeed,IncludeCellFunFilter, FunScore, yaml, confSet,CaData.statCell,NonTargetPath,NonTargets, IndexNonTargetTrial,CaData);
%%