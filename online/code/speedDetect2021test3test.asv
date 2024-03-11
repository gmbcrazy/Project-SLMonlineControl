%%
% --- Ana revision 08/22/2022 --- just check; nothing changed
% --- Ana revision 06/20/2022 --- tiny revision;
% --- Ana revision 06/16/2022 --- opts.DataRange = 'A2';
% --- Ana revision 06/13/2022 --- reverted to high precision w/o exclusions
% --- Ana revision 05/30/2022 --- no change; just check
% --- Ana revision 05/23/2022 --- removed excelSheetNum from input vars +
% some minor spellings + thresholds (exclusion of tiny events)
% --- Ana revision 05/12/2022 --- xlsread replacement
% --- Ana 11/19/2021 ---


function [fTrial, fTrialSpeed] = speedDetect2021test3test(processedDir, excelfilename, rowsOfInterest)
processedDir        = 'E:\Hari\11232022\statest\';%'E:\Hari\11242022\statest\';
excelfilename       = 'C:\test\2022_VIPSSTSummaryCOPY.xlsx';
rowsOfInterest      = [18];%18
%[fTrial, fTrialBehav1, fTrialSpeed] = speedDetect2021(processedDir,
%excelfilename, rowsOfInterest)%original
% --- Load vars ---
expFiles = dir(processedDir);
expFiles = expFiles(~strncmpi('.', {expFiles.name}, 1));

opts = detectImportOptions(excelfilename); opts.DataRange = 'A2';
text = readmatrix(excelfilename, opts);

header = opts.VariableNames;
for i = 1:numel(header)
    if ~isempty(header{i})
        header{i} = strrep(header{i},' ','');
        param.(header{i}) = text(1:end,i);
    end
end
shift = 1;
clear i opts header text

for m = 1:numel(rowsOfInterest)
    
    k = rowsOfInterest(m)-shift;
    
    disp(['---Processing Line' num2str(k+shift) '---'])
    
    currExp = expFiles(contains({expFiles.name}, param.Animal{k}) == 1 & contains({expFiles.name}, param.Date{k}) == 1 & contains({expFiles.name}, param.ExpDivision{k}) == 1);
    currExpPlane = currExp(1);
    currExpVars = dir([currExpPlane.folder '\' currExpPlane.name '\*vars.mat']);
    behavVars = dir([currExpPlane.folder '\' currExpPlane.name '\*Events.mat']);
    groupVars = dir([currExpPlane.folder '\' currExpPlane.name '\*groups.mat']);%%added by Hari
    clear currExp currExpPlane
    
    load([currExpVars.folder '\' currExpVars.name], 'fTrial', 'fR','iscell','stat','ops')
% load([behavVars.folder '\' behavVars.name], 'fTrialBehav1')%By Hari
%    load([currExpVars.folder '\' groupVars.name], 'groups')%% added by Hari
    
    clear behavVars
    
    %%
    fTrialSpeed = [];
    f = fields(fTrial);
     for t = 1:length(f)
        
        fieldname = char(f(t));
        
        % Calculate speed/LedCl individual bouts with maximum precision
        SpeedAllT(:,t)= fTrial.(fieldname).fSpeed-min(fTrial.(fieldname).fSpeed);
        FoFAllT=fTrial.(fieldname).fDeltaFoF;
        FoFAllTF(:,1:size(FoFAllT,2),t)=FoFAllT;
     end
        SpeedAllT=SpeedAllT(:);
        peaks = find(SpeedAllT > 0.1); %0.002% need accuracy to start with, else I won't have precise onsets and offsets
        dt = [0; diff(peaks)];
        shifted = dt < 15; % fR;%2
        shifted = [false;  shifted(2:end)];
        tags = diff([shifted; false]);
        starts = tags>0;
        ends = tags<0;
        fTrialSpeed.onsets1 = peaks(starts);
        fTrialSpeed.offsets1 = peaks(ends);
        %% By Hari to remove the small bout

        offset1_minus_onset1=[fTrialSpeed.offsets1]-[fTrialSpeed.onsets1];
        offset1_minus_onset1_idx=find(offset1_minus_onset1<30); %% remove the bout of one second time window
      
        onset1rev=fTrialSpeed.onsets1 ;
        offset1rev=fTrialSpeed.offsets1 ;
       onset1rev(offset1_minus_onset1_idx)=[];
        offset1rev(offset1_minus_onset1_idx)=[];

        fTrialSpeed.onsets1 = onset1rev;
        fTrialSpeed.offsets1 = offset1rev;
        % Join bouts that are less than 0.5s part ("behavioral bout")
%         peaks = find(temp > 0.008);%%0.002
%         dt = [0; diff(peaks)];
%         shifted = dt < round(fR/2);
%         shifted = [false;  shifted(2:end)];
%         
%         tags = diff([shifted; false]);
%         starts = tags>0;
%         ends = tags<0;
%         fTrialSpeed.onsets2 = peaks(starts);
%         fTrialSpeed.offsets2 = peaks(ends);
        clear peaks dt shifted tags starts ends
        
        
    
    [~, new, ~] = fileparts([currExpVars.folder '\' currExpVars.name]);
    k = strfind(new,'_vars');
    new = new(1:k);
    %filename=([currExpVars.folder '\' new 'behavEvents.mat']);
   % save([currExpVars.folder '\' new 'behavEvents.mat'], 'fTrialSpeed', '-append')   
   save([currExpVars.folder '\' new 'behavEvents.mat'], 'fTrialSpeed','SpeedAllT')   
end
figure;
plot(SpeedAllT)
hold on
xline(onset1rev, 'g');
xline(offset1rev, 'r');
%%%%%%%%%%%%%%%%%added to this code by Hari
% onW1=60;
% offW1=60;
%  onset1revF=find(onset1rev<onW1);
%  onset1rev(onset1revF)=[];
% offset1revF=find(offset1rev+offW1>size(SpeedAllT,1));
% offset1rev(offset1revF)=[];

siti=1;
itiW1=150;
for s=1:length(onset1rev)-1
    tempItI=onset1rev(s+1)-offset1rev(s);
    if tempItI>itiW1
    onset1reviti(siti)=onset1rev(s);
    offset1reviti(siti)=offset1rev(s);
    siti=siti+1;
    end
end
% figure;
% plot(SpeedAllT)
% hold on
% xline(onset1reviti, 'g');
% xline(offset1reviti, 'r');

onset1reviti=onset1rev(1:end);
offset1reviti=offset1rev(1:end);
%%%%%%%%%%%%%%%%%%%%%
FoFAllTFp=permute(FoFAllTF,[1 3 2]);
FoFAllTF= reshape(FoFAllTFp,[],size(FoFAllTF,2),1);
iscell = iscell; ce = find(iscell(:,1) == 1); 
FoFAllTF=FoFAllTF(:,ce);
FoFAllTF=FoFAllTF';
OnW1=30*2;%%30 =1 second therefore fr*1=1sec
OnW1plot=30*5;
OfW1=30*3;%%30 =1 second
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for event=1:length(onset1reviti)%cellN=1:size(FoFAllTF,1)
 for cellN=1:size(FoFAllTF,1)%event=1:length(onset1reviti)
     cellTriggertemp=FoFAllTF(cellN,onset1reviti(event)-OnW1:offset1reviti(event)-1);
     cellTriggertempplot=FoFAllTF(cellN,onset1reviti(event)-OnW1plot:offset1reviti(event)-1);

     cellTriggertempNorm(cellN,1,event)=[sum(cellTriggertemp)]/[numel(cellTriggertemp)/fR];%

%cellbasetemp=FoFAllTF(cellN,onset1reviti(event)-OnW1plot:onset1reviti(event)-OnW1plot+OnW1-1);
cellbasetemp=FoFAllTF(cellN,onset1reviti(event)-60:onset1reviti(event)-OnW1-1);

     cellbasetempNorm(cellN,1,event)=[sum(cellbasetemp)]/[numel(cellbasetemp)/fR];%
     %cellbasetemp=cellbasetemp(length(cellTriggertemp)+1:end);
     speedTrigplot=SpeedAllT(onset1reviti(event)-OnW1plot:offset1reviti(event)-1);
     speedTrigTplot(event,1:length(speedTrigplot))=speedTrigplot;
%%plot
     speedTrig=SpeedAllT(onset1reviti(event)-OnW1:offset1reviti(event)-1);
     speedTrigT(event,1:length(speedTrig))=speedTrig;
%%
     FoFAllTFc=FoFAllTF(cellN,:);
     cal_speed_corr(cellN)=corr(FoFAllTFc(:),SpeedAllT);

      FoFTrig(cellN,1:length(cellTriggertemp),event)=cellTriggertemp;
      FoFTrigplot(cellN,1:length(cellTriggertempplot),event)=cellTriggertempplot;
     FoFTrigbase(cellN,1:length(cellbasetemp),event)=cellbasetemp;
 end
end
%%plot
FoFTrigmean=nanmean(FoFTrig,3);
speedTrigTmean=nanmean(speedTrigT);

h = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(10,1,[1 6])
climits = [0 10];
imagesc(FoFTrigmean,climits)
c = colorbar;
pos = get(c,'Position');
c.Position = [pos(1)+.04, pos(2), 0.01, pos(4)]; % move colorbar so graphs are aligned
c.Label.String = 'dF scale)';
%colormap bone
title('\fontsize{11}dFoF')
ylabel('Cell #')
xticks(3600);
xticklabels('2 min')

subplot(10,1,[8 10])
timeScaleFOF=(-fR*2:size(FoFTrigmean,2))/fR;%%need to change here if we change the negative from onset i.e. fr*2 if negative from onset is 2
plot(timeScaleFOF(1:length(speedTrigTmean)),speedTrigTmean, 'Color', [.85 .33 .1], 'LineWidth', 1)
xlim([-2 max(timeScaleFOF)]);%% change here also
ylim([0 max(speedTrigTmean)]);
title('\fontsize{10}Locomotion')
ylabel('\fontsize{10}Speed')
xlabel('\fontsize{10}time (second)')
box off

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for t=1:size(cellTriggertempNorm,1)
    cellTriggertempNorm_cell1=cellTriggertempNorm(t,1,:);
    cellbasetempNorm_cell1=cellbasetempNorm(t,1,:);
    [~, Trigger_baseP(t)]=ttest(cellTriggertempNorm_cell1(:),cellbasetempNorm_cell1(:),"Alpha",0.05);
end
Trigger_baseP=Trigger_baseP(:);
[ cell_Provalue cell_idx]=sort(Trigger_baseP);
cell_Provalue_thr=cell_Provalue(find(cell_Provalue<0.05));
cell_idx_thr=cell_idx(1:numel(cell_Provalue_thr));
cell_R_thr=cal_speed_corr(cell_idx_thr);
T_max2min=table(cell_idx_thr(:),cell_Provalue_thr(:),cell_R_thr(:),'VariableNames',{'cellID','Pvalue','R'});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FoFTrigmean=nanmean(FoFTrigplot,3);
FoFTrigmean_statSing=FoFTrigmean(cell_idx_thr,:);
%%
[ cell_r2 r2_idx]=sort(cell_R_thr,'descend');
FoFTrigmean_statSingRbased=FoFTrigmean_statSing(r2_idx,:);
%%
speedTrigTmeanplot=nanmean(speedTrigTplot);

h = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(10,1,[1 6])
climits = [0 10];
imagesc(FoFTrigmean_statSing,climits)
c = colorbar;
pos = get(c,'Position');
c.Position = [pos(1)+.04, pos(2), 0.01, pos(4)]; % move colorbar so graphs are aligned
c.Label.String = 'dF scale)';
%colormap bone
title('\fontsize{11}dFoF')
ylabel('Cell #')
xticks(3600);
xticklabels('2 min')

subplot(10,1,[8 10])
%timeScaleFOF=(1:size(FoFTrigmean,2))/30;
timeScaleFOF2=(-fR*5:size(FoFTrigmean,2))/fR;
plot(timeScaleFOF2(1:length(speedTrigTmeanplot)),speedTrigTmeanplot, 'Color', [.85 .33 .1], 'LineWidth', 1)
xlim([min(timeScaleFOF2) max(timeScaleFOF2)]);
ylim([0 max(speedTrigTmeanplot)]);
title('\fontsize{10}Locomotion')
ylabel('\fontsize{10}Speed')
xlabel('\fontsize{10}time (second)')
box off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot(cell_R_thr, cell_Provalue_thr, 'b.', 'MarkerSize', 15);
lsline
ylabel('pvalue')
xlabel('R')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FoFAllTF
% cell_idx_thr
%
h1= figure('units','normalized','outerposition',[0 0 1 1]);
hold on
% subplot(10,1,[1 6])
%SpeedAllT
%%
[ cell_r2wp r2_idxwp]=sort(cal_speed_corr,'descend');%%without p
%%
xtime=1:1:size(FoFAllTF,2);
xtime=xtime/30;
for p=1:numel(cell_idx)%r2_idxwp
  cell_idx_thrtemp=cell_idx(p);
 subplot(numel(cell_idx),1,p)
 plot(xtime,FoFAllTF(cell_idx_thrtemp,:), 'Color', [.85 .33 .1], 'LineWidth', 1)
FoFAllTF_idx=FoFAllTF(cell_idx_thrtemp,:);
FoFAllTF_idx2=[FoFAllTF_idx-nanmean(FoFAllTF_idx)]/(std(FoFAllTF_idx));

end
hold on
subplot(numel(cell_idx)+1,1,p+1)
plot(xtime,SpeedAllT, 'Color', [.1 .1 .1], 'LineWidth', 1)
ylabel('\fontsize{8}Speed')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ci=1:length(stat)
    statiscell=stat{ci}.med;
   statiscellTotal(ci,:)=statiscell;
end
%cell_idx
%cal_speed_corrf=find(cal_speed_corr>0.5);
%ROIthr=statiscellTotal(cal_speed_corrf,:);
ROIthr=statiscellTotal(cell_idx_thr,:);
Plane_fileXYcorrTotal=ROIthr;
%%for Naparm points
% Plane_fileXYcorrTotal=flip(Plane_fileXYcorrTotal,2);
% Plane_fileXYcorrTotal(:,1)=512-Plane_fileXYcorrTotal(:,1);
% Plane_fileXYcorrTotal=round(Plane_fileXYcorrTotal);
points.h=[];
points.X=(Plane_fileXYcorrTotal(:,2))';%%% xy coordinates are opposites in python suite2p
points.Y=(Plane_fileXYcorrTotal(:,1))';
points.Z=ones(1,size(Plane_fileXYcorrTotal,1));
points.Zum=99.9*(ones(1,size(Plane_fileXYcorrTotal,1)));
points.OffsetX=ones(1,size(Plane_fileXYcorrTotal,1));
points.OffsetY=ones(1,size(Plane_fileXYcorrTotal,1));
points.Idx=1:1:size(Plane_fileXYcorrTotal,1);
points.Img=ones(1,size(Plane_fileXYcorrTotal,1));
points.Group=nan(1,size(Plane_fileXYcorrTotal,1));
points.GroupCentroidX=nan(1,size(Plane_fileXYcorrTotal,1));
points.GroupCentroidY=nan(1,size(Plane_fileXYcorrTotal,1));
points.Counter=size(Plane_fileXYcorrTotal,1);
points.Weight=ones(1,size(Plane_fileXYcorrTotal,1));
points.Selected=zeros(1,size(Plane_fileXYcorrTotal,1));
%%Here Hari generate random point in each group for testing purpose 
roi_size=numel(points.Group);
roi_idx=1:roi_size;
if roi_size<50
roi_sort=1:roi_size;
roi_rand=roi_sort(sort(randperm(roi_size,2)));
roi_rand_diff=diff(roi_rand);
roi_rand_diff(end+1)=roi_size-sum(roi_rand_diff);
sum(roi_rand_diff);

elseif  roi_size>50 && roi_size<100
    roi_sort=1:roi_size;
roi_rand=roi_sort(sort(randperm(roi_size,12)));
roi_rand_diff=diff(roi_rand);
roi_rand_diff(end+1)=roi_size-sum(roi_rand_diff);
sum(roi_rand_diff);

else 
    roi_sort=1:roi_size;
    roi_rand=roi_sort(sort(randperm(roi_size,20)));
    roi_rand_diff=diff(roi_rand);
    roi_rand_diff(end+1)=roi_size-sum(roi_rand_diff);
    sum(roi_rand_diff);
end

roi_rand_diffCumsum=(cumsum(roi_rand_diff))+1;
for g=1:numel(roi_rand_diff)
    roi_sizeTemp=roi_idx(roi_rand_diffCumsum(g)-roi_rand_diff(g):roi_rand_diffCumsum(g)-1); 
    roi_group(roi_sizeTemp)=g;
end
points.Group=roi_group;%% to generate the groups of different points
save([processedDir,'POINT.mat'],'points')
points={'X','Y','Z','Zum','OffsetX','OffsetY','Idx','Img','Group','GroupCentroidX','GroupCentroidY','Counter','Weight','Selected'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%mean img

