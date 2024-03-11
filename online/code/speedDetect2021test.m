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


function [fTrial, fTrialSpeed] = speedDetect2021test(processedDir, excelfilename, rowsOfInterest)
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
FoFAllTFp=permute(FoFAllTF,[1 3 2]);
FoFAllTF= reshape(FoFAllTFp,[],size(FoFAllTF,2),1);
iscell = iscell; ce = find(iscell(:,1) == 1); 
FoFAllTF=FoFAllTF(:,ce);
FoFAllTF=FoFAllTF';
OnW1=30;%%30 =1 second
OfW1=90;%%30 =1 second
for event=1:length(onset1reviti)
 for cellN=1:size(FoFAllTF,1)
     cellTriggertemp=FoFAllTF(cellN,onset1reviti(event)-OnW1:offset1reviti(event)-1);
     cellbasetemp=FoFAllTF(cellN,onset1reviti(event)-OnW1:offset1reviti(event)+OfW1-1);
     cellbasetemp=cellbasetemp(length(cellTriggertemp)+1:end);
     speedTrig=SpeedAllT(onset1reviti(event)-OnW1:offset1reviti(event)-1);
     if length(cellTriggertemp)>length(cellbasetemp)
      sizeeql=numel(cellTriggertemp)-numel(cellbasetemp);
      cellbasetempQ=[cellbasetemp(:);nan(sizeeql,1)];
       [~, Act_basettestP]=ttest(cellTriggertemp(:),cellbasetempQ,"Alpha",0.01);
     elseif length(cellbasetemp)>length(cellTriggertemp)
              sizeeql=numel(cellbasetemp)-numel(cellTriggertemp);
             cellbasetempQ=[cellTriggertemp(:);nan(sizeeql,1)];
              [~, Act_basettestP]=ttest(cellbasetemp(:),cellbasetempQ,"Alpha",0.01);
     end
     cellTriggertemp=cellTriggertemp(:);
     %[~, Act_basettestP]=ttest(cellTriggertemp,cellbasetempQ,"Alpha",0.01);%paired t test
     Act_basettestPT(cellN,1:3,event)=[onset1reviti(event) offset1reviti(event) Act_basettestP];

     FoFSpeedR=corr(cellTriggertemp,speedTrig);
     FoFSpeedRT(cellN,1,event)=FoFSpeedR;

     FoFTrig(cellN,1:length(cellTriggertemp),event)=cellTriggertemp;
     FoFTrigbase(cellN,1:length(cellbasetemp),event)=cellbasetemp;
     speedTrigT(1:length(speedTrig),event)=speedTrig;
     
 end
end
%%
Act_basettestPTpvalue=Act_basettestPT(:,3,:);
for p=1:size(Act_basettestPTpvalue)
    Act_basettestPTpvaluecell1=Act_basettestPTpvalue(p,:);
    Act_basettestPTpvaluecell1f=find(Act_basettestPTpvaluecell1<0.01);
    Act_basettestPTpvaluecell1pro(p)=numel(Act_basettestPTpvaluecell1f)/size(Act_basettestPTpvalue,3);

end
pro_thr=0.7;
Act_basettestPTpvaluecell1pro_thr=find(Act_basettestPTpvaluecell1pro>pro_thr);
cell_ID_pro_thr=Act_basettestPTpvaluecell1pro_thr(:);
cell_Pro=Act_basettestPTpvaluecell1pro(cell_ID_pro_thr);

for r=1:size(FoFSpeedRT,1)
    speedRcell=FoFSpeedRT(r,1,:);
    speedRcell=speedRcell(:);
    speedFRcell=find(speedRcell>0);
    speedPRcellTt=speedRcell(speedFRcell);%%Totaltrial
    speedPRcellTtmean(r)=nanmean(speedPRcellTt);
end
cell_ID_R=speedPRcellTtmean(cell_ID_pro_thr);
T=table(cell_ID_pro_thr(:),cell_Pro(:),cell_ID_R(:),'VariableNames',{'cellID','Probability','R'});

%%second methods with max to min sort
[ cell_Provalue cell_idx]=sort(Act_basettestPTpvaluecell1pro,'descend');
cell_Provalue_thr=cell_Provalue(find(cell_Provalue>pro_thr));
cell_idx_thr=cell_idx(1:numel(cell_Provalue_thr));
cell_R_thr=speedPRcellTtmean(cell_idx_thr);
T_max2min=table(cell_idx_thr(:),cell_Provalue_thr(:),cell_R_thr(:),'VariableNames',{'cellID','Probability','R'});
%%
FoFTrigmean=nanmean(FoFTrig,3);
speedTrigTmean=nanmean(speedTrigT,2);

%figure;
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
timeScaleFOF=(1:size(FoFTrigmean,2))/30;
plot(timeScaleFOF,speedTrigTmean, 'Color', [.85 .33 .1], 'LineWidth', 1)
%xlim([1 size(speedTrigTmean,1)]);
ylim([0 max(speedTrigTmean)]);
title('\fontsize{10}Locomotion')
ylabel('\fontsize{10}Speed')
xlabel('\fontsize{10}time (second)')
box off
%xticks(3600);
%xticklabels('\fontsize{10}2 min')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%

x

