
clear all
AllPath{1}='\\nimhlabstore1.nimh.nih.gov\UFNC\FNC3Z\RawRecording\SL0777-Ai203\10292024\';
AllPath{end+1}='\\nimhlabstore1.nimh.nih.gov\UFNC\FNC3Z\RawRecording\SL0777-Ai203\10302014\';
AllPath{end+1}='\\nimhlabstore1.nimh.nih.gov\UFNC\FNC3Z\RawRecording\SL0242-Ai203\10012024\';
AllPath{end+1}='\\nimhlabstore1.nimh.nih.gov\UFNC\FNC3Z\RawRecording\SL0242-Ai203\09302024\';
AllPath{end+1}='\\nimhlabstore1.nimh.nih.gov\UFNC\FNC3Z\RawRecording\SL0242-Ai203\09042024\';
AllPath{end+1}='\\nimhlabstore1.nimh.nih.gov\UFNC\FNC3Z\RawRecording\SL0242-Ai203\08292024\';
AllPath{end+1}='\\nimhlabstore1.nimh.nih.gov\UFNC\FNC3Z\RawRecording\SL0543-Ai203\09102024\';






% LoadPath='\\nimhlabstore1.nimh.nih.gov\UFNC\FNC3Z\RawRecording\SL0777-Ai203\10292024\';

fs=6.9;

for i=1:length(AllPath)
    LoadPath=AllPath{i};
   [rSpeed{i},rStim{i},nFrame(i),fSpeedOri{i}]=FunScore_Cal(LoadPath,fs);
end


DisX=[-0.3:0.05:0.4];
DisSpeed=[0:0.1:3];
TopN=13;

%%
figure;
for i=1:length(AllPath)
    subplot(2,length(AllPath),i)
    m=histPlotLU(rSpeed{i},DisX,[0.5 0.1 0.1],0.5);
    hold on;
    r=sort(rSpeed{i},'descend');
    r(TopN);
    plot([r(TopN),r(TopN)],[0 max(m)],'Color',[0.1 0.1 0.9]);
    text(r(TopN),5,['Top' num2str(TopN)],'Color',[0.1 0.1 0.9])


    Ymax1=max(m)*0.75; 
    strInfo1=['Frame#' num2str(nFrame(i)), ' ' showNum(nFrame(i)/6.9/60,1) 'min'];
    text(0.05,Ymax1,strInfo1)

    Ymax1=max(m)*0.9; 
    strInfo1=['Cell# ' num2str(length(rSpeed{i}))];
    text(0.05,Ymax1,strInfo1)

    set(gca,'xtick',[-0.5:0.1:0.5],'xlim',[min(DisX) max(DisX)])
    if i==4
    xlabel('Corr-Speed');
    end

    if i==1
    ylabel('Cell#');
    end
    subplot(2,length(AllPath),i+length(AllPath))
    speedData=fSpeedOri{i}(:,1);
    speedData(speedData<0.1)=[];
    m=histPlotLU(speedData,DisSpeed,[0.5 0.1 0.1],0.5);
    Ymax=max(m)*0.75;
    if i==4
    xlabel('Speed Amp');
    end
    
    if i==1
    ylabel('Sample#');
    end
end