%% Initialize, run this part before mannual correction of Suite2p processed data for saving time
clear all
ConfigFolder='E:\LuSLMOnlineTest\SL0838-Ai203\01242025\';
load('C:\Users\zhangl33\Projects\Project-SLMonlineControl\CodeAndPackages\subfun\Color\colorMapPN3.mat')

ConfigFile='CurrentSLMsetting.yml';
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
TopCellN=14;  
SavePathStimSpeed=[SavePath 'Top' num2str(TopCellN) 'SpeedStimEdgeExc\'];
mkdir(SavePathStimSpeed);

step3_SubStep4_FuncFilterMPs;


 % A=load('E:\LuSLMOnlineTest\SL0838-Ai203\01242025\SingleP\Top19SpeedStimEdgeExc\AllIncluded\SLMIncludedIndFromIscell.mat');
 % A=load('E:\LuSLMOnlineTest\SL0838-Ai203\01242025\SingleP\Top19SpeedStimEdgeExc\AllIncludedOrigin\SLMIncludedIndFromIscell.mat');

%% After Mannual Adjust Generated Functional Filtered MarkPoints in PV
PostGPLadjusted(SavePathStimSpeed,IncludeCellFunFilter, FunScore, yaml, confSet,CaData.statCell,NonTargetPath,NonTargets, IndexNonTargetTrial,CaData);
%%
% Pos3D(A.SLMIncludedIndFromIscell,:)-A.Pos3Dneed



figure;
plot(rSpeed(CenterCloseI),'r-');
hold on;
plot(rStim(CenterCloseI),'g-');
plot(rCenterISpeed(1:TopCellN),0.3,'ro')
plot(rCenterIStim(1:TopCellN),0.3,'g*')



figure;
ImgDataTemp=[];
for iCell=1:length(IncludeCellFunFilter)
% subplot(TopCellN,1,i)
tempData=NeuroData(:,iscell(IncludeCellFunFilter(iCell)));
ImgDataTemp(iCell,:)=SmoothDec(tempData,0.1)';
end
imagesc(ImgDataTemp)
hold on;
plot(-mean(fSpeed,2)*(TopCellN+1),'r');
set(gca,'ylim',[-TopCellN TopCellN])
plot(-mean(fStim,2)*(TopCellN+1),'b');
set(gca,'ylim',[-TopCellN TopCellN])


figure;
tempData=[];
for iTop=1:length(rCenterISpeed)
% subplot(TopCellN,1,i)
tempData1=zscore(NeuroDataCell(:,rCenterISpeed(iTop)));
tempData(iTop,:)=SmoothDec(tempData1,3);
% plot(tempData+iTop);
end
fSpeedSmooth=SmoothDec(mean(fSpeed,2),1);
imagesc(tempData)
hold on;
colormap(gray)
set(gca,'clim',[-1 3])
plot(fSpeedSmooth*(length(rCenterISpeed)+1),'r');


figure;
for iTop=1:length(rCenterIStim)
% subplot(TopCellN,1,i)
tempData=NeuroDataCell(:,rCenterIStim(iTop));
tempData=SmoothDec(tempData,1);
plot(tempData+iTop);
hold on;
end
plot(mean(fStim,2)*(length(rCenterIStim)+1),'r');


subplot(2,1,1)
imagesc(NeuroData(:,iscell(CenterCloseI))');
subplot(2,1,2)
plot(mean(fSpeed,2));hold on;


figure;
for i=1:length(CenterCloseI)
    subplot(length(CenterCloseI),1,i)
plot(NeuroData(:,iscell(CenterCloseI(i))));hold on;
plot(mean(fStim,2)+5,'r-');
end


figure;
plot(rSpeed,'r');
hold on;
plot(rStim,'g');
hold on;
plot(rFspks,'b')

figure;
plot(-log10(pFspks),'b')
[~,r1]=sort(rFspks(CenterCloseI),'descend')  %%Cell locates close to edge of the view, were not considered as SML targets.
TopCellN=8;         %%   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<Top X cells highly associated with speed 
TopSpeedCellI=CenterCloseI(r1(1:TopCellN));
SavePathSpeed=[SavePath 'Top' num2str(TopCellN) 'SpeedEdgeExc\']
mkdir(SavePathSpeed)
XYZtoMarkPoint_NT(SavePathSpeed,Pos3D,TopSpeedCellI,NonTargets, IndexNonTargetTrial, yaml,confSet)





%% 
% load('C:\Users\zhangl33\Projects\GenMatCode\Plotfun\Color\colorMapPN3.mat');
% rCellspks=partialcorr(CaData.spks(iscell,:)',fSpeed);
rCellspks=corr(double(CaData.spks(iscell,:))','type','Spearman','rows','pairwise');
load('C:\Users\User\Project-SLMonlineControl\subfun\Color\colorMapPN3.mat')
AdjMat=rCellspks+1;  %%%%%%%%%%Non-negative weighted correlation.
clear SampleCorr;
for itt=1:size(AdjMat,1)
    AdjMat(itt,itt)=0;%%%%%%diagonal zeros for Adjcent matrix
end
k = full(sum(AdjMat));
twom = sum(k);
B = full(AdjMat - gamma*k'*k/twom);
tic
k = full(sum(AdjMat));
twom = sum(k);
gamma = 1.05;  %%Community Clustering Paramter%%
limit =100000; %%memory consideration for community clustering 
B = @(i) AdjMat(:,i) - gamma*k'*k(i)/twom;
disp('Clustering ...iterated_genlouvain.m');
%%%%%%Clustering
[SS,QQ,n_it]=iterated_genlouvain(B,limit,0);
max(SS)
toc
[~,r2] = sort(SS);
%         QQ = QQ/twom;
figure;
% imagesc(rCellspks(r2,r2))
AdjComImagesc(rCellspks,SS);
colormap(ColorPN3)
clim([-0.3 0.3])

IndSub=find(SS==3)
CellSub=double(CaData.spks(iscell(IndSub),:)');

[coeff,score,latent,tsquared] = pca(CellSub);
corr(score(:,1:3),fSpeed)



[W,H]=nnmf(double(CaData.spks(iscell,:))',3)
biplot(H','Scores',W);
axis([0 1.1 0 1.1])
xlabel('Column 1')
ylabel('Column 2')

r=corr(W,fSpeed)


% 1 2 3 7 9 13

clear Group
Group(1).Indices=[1:length(tempI)];
MarkPoints3D_GPLmaker(Pos3Dneed, yaml, true, SpiralSizeUM, SpiralRevolutions, SaveName,Group)
MarkPoints3D_XMLmaker_Points(Pos3Dneed,yaml,true, Repetition, SpiralSizeUM, SpiralRevolutions,UncagingLaserPower, SavePathAllPoint)


% MarkPoints3D_XMLmaker_Group(Group,Repetition, SpiralSizeUM, UncagingLaserPower, SavePath)


%  2 4 7 13 14 15 18 23 27

% tempI=[2 3 4 5 9 11 14 16 21 23]
tempI=[2 3 5 9 14 16 21 23]


