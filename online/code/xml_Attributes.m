function xml_Attributes()
%file=('C:\Work\2pdata\h04202022\BrightnessOverTime-04202022-1007-831\BrightnessOverTime-04202022-1007-831.xml');
file=('C:\Work\2pdata\h04222022\BrightnessOverTime-04222022-0951-034\BrightnessOverTime-04222022-0951-034.xml');
[s]=xml2struct( file );
Frame_time_series=s.PVScan.Sequence.Frame;
for i=1:size(Frame_time_series,2)
    Frame_time_series_absTime(i)=str2num(Frame_time_series{1,i}.Attributes.absoluteTime);
    %Frame_time_series_absTime_total(i)=Frame_time_series_absTime_temp;

    Frame_time_series_relTime(i)=str2num(Frame_time_series{1,i}.Attributes.relativeTime);
    %Frame_time_series_relTime_total(i)=Frame_time_series_relTime_temp;
end
Frame_time_series_abs_relTdiff=Frame_time_series_absTime-Frame_time_series_relTime;
%%%diff between abs and relative always same (first point,i.e. 1.6)
%figure;
plot(Frame_time_series_absTime,ones(size(Frame_time_series_absTime,2),1),'|b')
%hold on 
%plot(Frame_time_series_abs_relTdiff,ones(size(Frame_time_series_absTime,1),2),'|','r')
%frame_time_series=s.Children.Name;

%%% mark point xml to str (initial delay)

file=('C:\Work\2pdata\BrightnessOverTime-04202022-1007-830\BrightnessOverTime-04202022-1007-830_Cycle00001_MarkPoints.xml');
[s2]=xml2struct( file );
Ini_delay=str2num(s2.PVMarkPointSeriesElements.PVMarkPointElement.PVGalvoPointElement.Attributes.InitialDelay);
Ini_shift=(Ini_delay/1000)+Frame_time_series_absTime(1);%%change to seconds
%%%% read intensity over time exel file
[ITimenum,ITimetxt,ITimeraw] = xlsread('C:\Work\2pdata\h04222022\BrightnessOverTime-04222022-0951-034\BrightnessOverTime-04222022-0951-034_Cycle00001-botData') ;

figure;
plot(ITimenum(:,1),ITimenum(:,2),'r')
hold on
plot(ITimenum(:,1),ITimenum(:,3),'g')
hold on
plot(ITimenum(:,1),ITimenum(:,4),'b')
hold on
plot(Ini_shift,1:1200,'|k')




% filename=('C:\Work\2pdata\h04202022\BrightnessOverTime-04202022-1007-831\BrightnessOverTime-04202022-1007-831.xml');
% %filename=('C:\Work\2pdata\h04222022\BrightnessOverTime-04222022-0951-034\BrightnessOverTime-04222022-0951-034_Cycle00001_VoltageOutput_001.xml');
% theStruct = parseXML(filename);
x