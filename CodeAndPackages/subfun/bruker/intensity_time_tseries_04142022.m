function intensity_time_tseries_04142022(session)
%Dir_path=['C:\Work\2pdata\',session];
Dir_path=['D:\',session(1:end-2),'\',session(1:end)];
cd(Dir_path)
files = dir('T*.');
files_name=files.name;
cd(files_name)
files_exel=dir('T*.csv');

ITimenumT_total=nan(210,20);
ITimenumI_total=nan(210,20);
ic=1;
kc=1;
for k=1:2:size(files_exel,1)
file_name=files_exel(k).name;
ITimenum_odd = readmatrix(file_name);
ITimenum_odd_time(:,kc)=ITimenum_odd(1:59,1);
ITimenum_odd_I(:,kc)=ITimenum_odd(1:59,2);
kc=kc+1;
end
kc2=1;
for kk=2:2:size(files_exel,1)
file_name=files_exel(kk).name;
ITimenum_even = readmatrix(file_name);
ITimenum_even_time(:,kc2)=ITimenum_even(1:209,1);
ITimenum_even_I(:,kc2)=ITimenum_even(1:209,2);
kc2=kc2+1;
end

for kkk=1:size(ITimenum_odd_I,2)
    ITimenum_odd_even_Ibs=ITimenum_odd_I(:,kkk);
   ITimenum_odd_even_Ial=ITimenum_even_I(:,kkk);
ITimenum_odd_even_temp=[ITimenum_odd_even_Ibs;ITimenum_odd_even_Ial];
ITimenum_odd_even_total(:,kkk)=ITimenum_odd_even_temp;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ITimenum_mat_BS_cell2mat_time_trail_1 = 0:1/30:(268-1)/30; %t = 0:1/fs:(length(x)-1)/fs;%%converting fr to time series
figure;

% test=importdata('D:\h04142022\h0414202214\TSeries-04142022-0935-014\test.txt');
% ITimenum_odd_even_total=ITimenum_odd_even_total*3;
% testf=find(ITimenum_odd_even_total(:,1));
% ITimenum_odd_even_total(testf)=test;

ITimenum_mat_BS_cell2mat_I_trail_avg=nanmean(ITimenum_odd_even_total,2);
plot(ITimenum_mat_BS_cell2mat_time_trail_1,ITimenum_mat_BS_cell2mat_I_trail_avg,'k')
hold on
plot(2,1:max(ITimenum_mat_BS_cell2mat_I_trail_avg),'|k')

 hold on
hline = plot(ITimenum_mat_BS_cell2mat_time_trail_1,ITimenum_odd_even_total);
for j=1:length(hline)
    hline(j).Color = [hline(j).Color 0.4];  % alpha=0.1
end
hold on
xlabel('Time in seconds')
 hold on
 ylabel('Absolute intensity')
x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%










