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


function [fTrial, fTrialSpeed] = speedDetect2021test2(processedDir, excelfilename, rowsOfInterest)
processedDir        = 'C:\Data\Hari\statest3\';
excelfilename       = 'C:\test\2022_VIPSSTSummaryCOPY.xlsx';
rowsOfInterest      = [17];
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
    
    load([currExpVars.folder '\' currExpVars.name], 'fTrial', 'fR','iscell')
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
        peaks = find(SpeedAllT > 0.5); %0.002% need accuracy to start with, else I won't have precise onsets and offsets
        dt = [0; diff(peaks)];
        shifted = dt < 2; % fR;%2
        shifted = [false;  shifted(2:end)];% shifted = [false;  shifted(2:end)];
        tags = diff([shifted; false]);
        starts = tags>0;
        ends = tags<0;
        fTrialSpeed.onsets1 = peaks(starts);
        fTrialSpeed.offsets1 = peaks(ends);
        % Join bouts that are less than 0.5s part ("behavioral bout")
        peaks = find(SpeedAllT > 0.1);
        dt = [0; diff(peaks)];
        shifted = dt < round(fR/2);
        shifted = [false;  shifted(2:end)];
        tags = diff([shifted; false]);
        starts = tags>0;
        ends = tags<0;
        fTrialSpeed.(fieldname).onsets2 = peaks(starts);
        fTrialSpeed.(fieldname).offsets2 = peaks(ends);
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
figure;
plot(SpeedAllT)
hold on
xline(onset1reviti, 'g');
xline(offset1reviti, 'r');
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

    
    cellbasetemp=FoFAllTF(cellN,onset1reviti(event)-OnW1-OfW1:onset1reviti(event)-OnW1-1);

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for t=1:size(cellTriggertempNorm,1)
    cellTriggertempNorm_cell1=cellTriggertempNorm(t,1,:);
    cellbasetempNorm_cell1=cellbasetempNorm(t,1,:);
    [~, Trigger_baseP(t)]=ttest(cellTriggertempNorm_cell1(:),cellbasetempNorm_cell1(:),"Alpha",0.01);
end
Trigger_baseP=Trigger_baseP(:);
[ cell_Provalue cell_idx]=sort(Trigger_baseP);
cell_Provalue_thr=cell_Provalue(find(cell_Provalue<0.05));
cell_idx_thr=cell_idx(1:numel(cell_Provalue_thr));
cell_R_thr=cal_speed_corr(cell_idx_thr);
T_max2min=table(cell_idx_thr(:),cell_Provalue_thr(:),cell_R_thr(:),'VariableNames',{'cellID','Pvalue','R'});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FoFTrigmean=nanmean(FoFTrigplot,3);
FoFTrigmean_statSing=FoFTrigmean(cell_idx_thr,:);
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
for p=1:numel(cell_idx_thr)
  cell_idx_thrtemp=cell_idx_thr(p);
 subplot(numel(cell_idx_thr),1,p)
 plot(FoFAllTF(cell_idx_thrtemp,:), 'Color', [.85 .33 .1], 'LineWidth', 1)
FoFAllTF_idx=FoFAllTF(cell_idx_thrtemp,:);
FoFAllTF_idx2=[FoFAllTF_idx-nanmean(FoFAllTF_idx)]/(std(FoFAllTF_idx));

end
hold on
subplot(numel(cell_idx_thr)+1,1,p+1)
plot(SpeedAllT, 'Color', [.1 .1 .1], 'LineWidth', 1)
ylabel('\fontsize{8}Speed')
%%
x