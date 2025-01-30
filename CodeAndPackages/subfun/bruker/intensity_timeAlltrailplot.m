function intensity_timeAlltrailplot(session)
Dir_path=['C:\Work\2pdata\',session];
%Dir_path=['D:\',session(1:end-2),'\',session(1:end)];
cd(Dir_path)
files = dir('B*.');

%ITimenum=nan(10000,10);

for i=1:size(files,1)
    file_name=files(i).name;
    cd(file_name)
%%intensity over time exl file
ITimenum = readmatrix([file_name,'_Cycle00001-botData','.csv']);
ITimenum_mat(i)={ITimenum};

%%
xml2strconver=xml2struct( [file_name,'.xml'] );
Frame_time_series=xml2strconver.PVScan.Sequence.Frame;

xml2strconverMarkP=xml2struct( [file_name,'_Cycle00001_MarkPoints.xml'] );%%All mark points has same delay i.e. 5000
MPIni_delay=str2num(xml2strconverMarkP.PVMarkPointSeriesElements.PVMarkPointElement.PVGalvoPointElement.Attributes.InitialDelay);
MPIni_delay=MPIni_delay/1000; %%mSecond to second
%%
    for j=1:size(Frame_time_series,2)
    Frame_time_series_absTime(i,j)=str2num(Frame_time_series{1,j}.Attributes.absoluteTime);
    Frame_time_series_relTime(i,j)=str2num(Frame_time_series{1,j}.Attributes.relativeTime);
    
    end
   cd(Dir_path)
end
ITimenum_mat_maxLength=max( cellfun(@length, ITimenum_mat));
ITimenum_mat_eqlSize=cellfun( @(v) padarray(v,[ITimenum_mat_maxLength-length(v),0],nan,'post'), ITimenum_mat, 'uni',false );
ITimenum_mat_cell2mat=cell2mat(ITimenum_mat_eqlSize);



% for ii=1:size(ITimenum_mat_eqlSize,2)
%     ITimenum_mat_cell2mat_idx=ITimenum_mat_eqlSize{ii};
%     ITimenum_mat_cell2mat_time(:,ii)=ITimenum_mat_cell2mat_idx(:,1);
%     ITimenum_mat_cell2mat_Int(:,ii)=ITimenum_mat_cell2mat_idx(:,2);
% end

if size(ITimenum_mat_eqlSize{1},2) == 2
ITimenum_mat_cell2mat_time=ITimenum_mat_cell2mat(:,1:2:end);
ITimenum_mat_cell2mat_time_1=ITimenum_mat_cell2mat_time(1,:);
ITimenum_mat_cell2mat_time_ral=ITimenum_mat_cell2mat_time-ITimenum_mat_cell2mat_time_1;
ITimenum_mat_cell2mat_Int=ITimenum_mat_cell2mat(:,2:2:end);

ITimenum_mat_cell2mat_time_avg=nanmean(ITimenum_mat_cell2mat_time_ral,2);
ITimenum_mat_cell2mat_Int_avg=nanmean(ITimenum_mat_cell2mat_Int,2);
figure;
plot(ITimenum_mat_cell2mat_time_avg,ITimenum_mat_cell2mat_Int_avg,'k')
hold on
plot(MPIni_delay,1:max(ITimenum_mat_cell2mat_Int_avg),'|k')
hold on
hline1 = plot(ITimenum_mat_cell2mat_time_avg,ITimenum_mat_cell2mat_Int);
for ij11=1:length(hline1)
    hline1(ij11).Color = [hline1(ij11).Color 0.4];  % alpha=0.3
end

% hold on
% xlabel('Time in Seconds')
% hold on 
% ylabel('Intensity')
 elseif size(ITimenum_mat_eqlSize{1},2) == 3
ITimenum_mat_cell2mat_time=ITimenum_mat_cell2mat(:,1:3:end);
ITimenum_mat_cell2mat_time_1=ITimenum_mat_cell2mat_time(1,:);
ITimenum_mat_cell2mat_time_ral=ITimenum_mat_cell2mat_time-ITimenum_mat_cell2mat_time_1;
ITimenum_mat_cell2mat_Int23=ITimenum_mat_cell2mat(:,2:3:end);
ITimenum_mat_cell2mat_Int33=ITimenum_mat_cell2mat(:,3:3:end);

ITimenum_mat_cell2mat_time_avg=nanmean(ITimenum_mat_cell2mat_time_ral,2);
ITimenum_mat_cell2mat_Int_avg23=nanmean(ITimenum_mat_cell2mat_Int23,2);
ITimenum_mat_cell2mat_Int_avg33=nanmean(ITimenum_mat_cell2mat_Int33,2);
figure;
subplot(2,1,1)
plot(ITimenum_mat_cell2mat_time_avg,ITimenum_mat_cell2mat_Int_avg23)
hold on
plot(MPIni_delay,1:max(ITimenum_mat_cell2mat_Int_avg23),'|k')
subplot(2,1,2)
plot(ITimenum_mat_cell2mat_time_avg,ITimenum_mat_cell2mat_Int_avg33)
hold on
plot(MPIni_delay,1:max(ITimenum_mat_cell2mat_Int_avg33),'|k')
% xlabel('Time in Seconds')
% hold on 
% ylabel('Mean Intensity')
elseif size(ITimenum_mat_eqlSize{1},2) == 4
ITimenum_mat_cell2mat_time=ITimenum_mat_cell2mat(:,1:4:end);
ITimenum_mat_cell2mat_time_1=ITimenum_mat_cell2mat_time(1,:);
ITimenum_mat_cell2mat_time_ral=ITimenum_mat_cell2mat_time-ITimenum_mat_cell2mat_time_1;
ITimenum_mat_cell2mat_Int22=ITimenum_mat_cell2mat(:,2:4:end);
ITimenum_mat_cell2mat_Int33=ITimenum_mat_cell2mat(:,3:4:end);
ITimenum_mat_cell2mat_Int44=ITimenum_mat_cell2mat(:,4:4:end);

ITimenum_mat_cell2mat_time_avg=nanmean(ITimenum_mat_cell2mat_time_ral,2);
ITimenum_mat_cell2mat_Int_avg22=nanmean(ITimenum_mat_cell2mat_Int22,2);
ITimenum_mat_cell2mat_Int_avg33=nanmean(ITimenum_mat_cell2mat_Int33,2);
ITimenum_mat_cell2mat_Int_avg44=nanmean(ITimenum_mat_cell2mat_Int44,2);
%%%%%%%%%%%%%%peak intensity
ITimenum_mat_cell2mat_Int22_peak=ITimenum_mat_cell2mat_Int22(160:210,:);
ITimenum_mat_cell2mat_Int22_peakMM=max(ITimenum_mat_cell2mat_Int22_peak)-min(ITimenum_mat_cell2mat_Int22_peak);
ITimenum_mat_cell2mat_Int22_peakMMmean=mean(ITimenum_mat_cell2mat_Int22_peakMM)
ITimenum_mat_cell2mat_Int22_peakMMstd=std(ITimenum_mat_cell2mat_Int22_peakMM)

ITimenum_mat_cell2mat_Int33_peak=ITimenum_mat_cell2mat_Int33(160:210,:);
ITimenum_mat_cell2mat_Int33_peakMM=max(ITimenum_mat_cell2mat_Int33_peak)-min(ITimenum_mat_cell2mat_Int33_peak);
ITimenum_mat_cell2mat_Int33_peakMMmean=mean(ITimenum_mat_cell2mat_Int33_peakMM)
ITimenum_mat_cell2mat_Int33_peakMMstd=std(ITimenum_mat_cell2mat_Int33_peakMM)

ITimenum_mat_cell2mat_Int44_peak=ITimenum_mat_cell2mat_Int44(160:210,:);
ITimenum_mat_cell2mat_Int44_peakMM=max(ITimenum_mat_cell2mat_Int44_peak)-min(ITimenum_mat_cell2mat_Int44_peak);
ITimenum_mat_cell2mat_Int44_peakMMmean=mean(ITimenum_mat_cell2mat_Int44_peakMM)
ITimenum_mat_cell2mat_Int44_peakMMstd=std(ITimenum_mat_cell2mat_Int44_peakMM)
%%%%%%%%%%%%%%
ITimenum_mat_cell2mat_Int22h=ITimenum_mat_cell2mat_Int22(1:300,:);
ITimenum_mat_cell2mat_Int33h=ITimenum_mat_cell2mat_Int33(1:300,:);
ITimenum_mat_cell2mat_Int44h=ITimenum_mat_cell2mat_Int44(1:300,:);

figure;
subplot(3,1,1)
hline22 = plot(ITimenum_mat_cell2mat_time_avg(1:300),ITimenum_mat_cell2mat_Int22h);
hold on 
plot(ITimenum_mat_cell2mat_time_avg(1:300),nanmean(ITimenum_mat_cell2mat_Int22h,2),'k')
hold on
plot(MPIni_delay,1:max(ITimenum_mat_cell2mat_Int_avg22),'|k')
hold on
ylabel('Absolute intensity')
%hline33= plot(ITimenum_mat_cell2mat_time_avg(1:300),ITimenum_mat_cell2mat_Int33h);
%hline44 = plot(ITimenum_mat_cell2mat_time_avg(1:300),ITimenum_mat_cell2mat_Int44h);
for ij22=1:length(hline22)
    hline22(ij22).Color = [hline22(ij22).Color 0.4];  % alpha=0.3
end
subplot(3,1,2)
hline33= plot(ITimenum_mat_cell2mat_time_avg(1:300),ITimenum_mat_cell2mat_Int33h);
hold on 
plot(ITimenum_mat_cell2mat_time_avg(1:300),nanmean(ITimenum_mat_cell2mat_Int33h,2),'k')
hold on
plot(MPIni_delay,1:max(ITimenum_mat_cell2mat_Int_avg33),'|k')
hold on
ylabel('Absolute intensity')
for ij33=1:length(hline33)
    hline33(ij33).Color = [hline33(ij33).Color 0.4];  % alpha=0.3
end

subplot(3,1,3)
hline44= plot(ITimenum_mat_cell2mat_time_avg(1:300),ITimenum_mat_cell2mat_Int44h);
hold on 
plot(ITimenum_mat_cell2mat_time_avg(1:300),nanmean(ITimenum_mat_cell2mat_Int44h,2),'k')
hold on
plot(MPIni_delay,1:max(ITimenum_mat_cell2mat_Int_avg44),'|k')
hold on
ylabel('Absolute intensity')
hold on
xlabel('Time in Seconds')
for ij44=1:length(hline44)
    hline44(ij44).Color = [hline33(ij44).Color 0.4];  % alpha=0.3
end


figure;
subplot(3,1,1)
plot(ITimenum_mat_cell2mat_time_avg(1:300),ITimenum_mat_cell2mat_Int_avg22(1:300))
hold on
plot(MPIni_delay,1:max(ITimenum_mat_cell2mat_Int_avg22),'|k')
hold on
ylabel('Absolute intensity')
subplot(3,1,2)
plot(ITimenum_mat_cell2mat_time_avg(1:300),ITimenum_mat_cell2mat_Int_avg33(1:300))
hold on
plot(MPIni_delay,1:max(ITimenum_mat_cell2mat_Int_avg33),'|k')
hold on
ylabel('Absolute intensity')
subplot(3,1,3)
plot(ITimenum_mat_cell2mat_time_avg(1:300),ITimenum_mat_cell2mat_Int_avg44(1:300))
hold on
plot(MPIni_delay,1:max(ITimenum_mat_cell2mat_Int_avg44),'|k')
hold on
ylabel('Absolute intensity')
hold on
xlabel('Time in Seconds')
% xlabel('Time in Seconds')
% hold on 
% ylabel('Mean Intensity')
elseif size(ITimenum_mat_eqlSize{1},2) == 5
ITimenum_mat_cell2mat_time=ITimenum_mat_cell2mat(:,1:5:end);
ITimenum_mat_cell2mat_time_1=ITimenum_mat_cell2mat_time(1,:);
ITimenum_mat_cell2mat_time_ral=ITimenum_mat_cell2mat_time-ITimenum_mat_cell2mat_time_1;
ITimenum_mat_cell2mat_Int22=ITimenum_mat_cell2mat(:,2:5:end);
ITimenum_mat_cell2mat_Int33=ITimenum_mat_cell2mat(:,3:5:end);
ITimenum_mat_cell2mat_Int44=ITimenum_mat_cell2mat(:,4:5:end);
ITimenum_mat_cell2mat_Int55=ITimenum_mat_cell2mat(:,5:5:end);



ITimenum_mat_cell2mat_time_avg=nanmean(ITimenum_mat_cell2mat_time_ral,2);
ITimenum_mat_cell2mat_Int_avg22=nanmean(ITimenum_mat_cell2mat_Int22,2);
ITimenum_mat_cell2mat_Int_avg33=nanmean(ITimenum_mat_cell2mat_Int33,2);
ITimenum_mat_cell2mat_Int_avg44=nanmean(ITimenum_mat_cell2mat_Int44,2);
ITimenum_mat_cell2mat_Int_avg55=nanmean(ITimenum_mat_cell2mat_Int55,2);

figure;
subplot(4,1,1)
plot(ITimenum_mat_cell2mat_time_avg,ITimenum_mat_cell2mat_Int_avg22)
hold on
plot(MPIni_delay,1:max(ITimenum_mat_cell2mat_Int_avg22),'|k')
subplot(4,1,2)
plot(ITimenum_mat_cell2mat_time_avg,ITimenum_mat_cell2mat_Int_avg33)
hold on
plot(MPIni_delay,1:max(ITimenum_mat_cell2mat_Int_avg33),'|k')
subplot(4,1,3)
plot(ITimenum_mat_cell2mat_time_avg,ITimenum_mat_cell2mat_Int_avg44)
hold on
plot(MPIni_delay,1:max(ITimenum_mat_cell2mat_Int_avg44),'|k')
subplot(4,1,4)
plot(ITimenum_mat_cell2mat_time_avg,ITimenum_mat_cell2mat_Int_avg55)
% hold on
% plot(MPIni_delay,1:max(ITimenum_mat_cell2mat_Int_avg55),'|k')
% xlabel('Time in Seconds')
% hold on 
% ylabel('Mean Intensity')
elseif size(ITimenum_mat_eqlSize{1},2) == 6
ITimenum_mat_cell2mat_time=ITimenum_mat_cell2mat(:,1:5:end);
ITimenum_mat_cell2mat_time_1=ITimenum_mat_cell2mat_time(1,:);
ITimenum_mat_cell2mat_time_ral=ITimenum_mat_cell2mat_time-ITimenum_mat_cell2mat_time_1;
ITimenum_mat_cell2mat_Int22=ITimenum_mat_cell2mat(:,2:5:end);
ITimenum_mat_cell2mat_Int33=ITimenum_mat_cell2mat(:,3:5:end);
ITimenum_mat_cell2mat_Int44=ITimenum_mat_cell2mat(:,4:5:end);
ITimenum_mat_cell2mat_Int55=ITimenum_mat_cell2mat(:,5:5:end);
ITimenum_mat_cell2mat_Int66=ITimenum_mat_cell2mat(:,6:5:end);


ITimenum_mat_cell2mat_time_avg=nanmean(ITimenum_mat_cell2mat_time_ral,2);
ITimenum_mat_cell2mat_Int_avg22=nanmean(ITimenum_mat_cell2mat_Int22,2);
ITimenum_mat_cell2mat_Int_avg33=nanmean(ITimenum_mat_cell2mat_Int33,2);
ITimenum_mat_cell2mat_Int_avg44=nanmean(ITimenum_mat_cell2mat_Int44,2);
ITimenum_mat_cell2mat_Int_avg55=nanmean(ITimenum_mat_cell2mat_Int55,2);
ITimenum_mat_cell2mat_Int_avg66=nanmean(ITimenum_mat_cell2mat_Int66,2);

figure;
subplot(5,1,1)
plot(ITimenum_mat_cell2mat_time_avg,ITimenum_mat_cell2mat_Int_avg22)
hold on
plot(MPIni_delay,1:max(ITimenum_mat_cell2mat_Int_avg22),'|k')
subplot(5,1,2)
plot(ITimenum_mat_cell2mat_time_avg,ITimenum_mat_cell2mat_Int_avg33)
hold on
plot(MPIni_delay,1:max(ITimenum_mat_cell2mat_Int_avg33),'|k')
subplot(5,1,3)
plot(ITimenum_mat_cell2mat_time_avg,ITimenum_mat_cell2mat_Int_avg44)
hold on
plot(MPIni_delay,1:max(ITimenum_mat_cell2mat_Int_avg44),'|k')
subplot(5,1,4)
plot(ITimenum_mat_cell2mat_time_avg,ITimenum_mat_cell2mat_Int_avg55)
hold on
plot(MPIni_delay,1:max(ITimenum_mat_cell2mat_Int_avg55),'|k')
subplot(5,1,5)
plot(ITimenum_mat_cell2mat_time_avg,ITimenum_mat_cell2mat_Int_avg66)
hold on
plot(MPIni_delay,1:max(ITimenum_mat_cell2mat_Int_avg66),'|k')
xlabel('Time in Seconds')
hold on 
ylabel('Mean Intensity')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%











x