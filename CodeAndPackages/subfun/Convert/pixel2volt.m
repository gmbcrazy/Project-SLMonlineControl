function [str3_xTotal str3_yTotal]=pixel2volt()
%voltage_x=importdata('D:\test512_512X.gpl');
voltage_x=importdata('D:\07222022\x_512.gpl');
str1_x = regexprep(voltage_x,'[,;=]', ' ');
str2_x = regexprep(regexprep(str1_x,'[^- 0-9.eE(,)]',''), ' \D* ',' ');
str3_x = regexprep(str2_x, {'\.\s','\E\s','\e\s','\s\E','\s\e'},' ');

%voltage_y=importdata('D:\test512_512Y.gpl');
voltage_y=importdata('D:\07222022\y_512.gpl');
str1_y= regexprep(voltage_y,'[,;=]', ' ');
str2_y = regexprep(regexprep(str1_y,'[^- 0-9.eE(,)]',''), ' \D* ',' ');
str3_y = regexprep(str2_y, {'\.\s','\E\s','\e\s','\s\E','\s\e'},' ');

%numArray = str2num(str3_x);
xyxount=1;
for xy=3:1:514%%1027
    str3_xtemp=str2num(str3_x{xy});
    str3_xtempTotal(xyxount)=str3_xtemp(1);
 
    str3_ytemp=str2num(str3_y{xy});
    str3_ytempTotal(xyxount)=str3_ytemp(2);
xyxount=xyxount+1;
end
str3_xTotal=str3_xtempTotal*2;
str3_yTotal=str3_ytempTotal*2;%%(3:2:end-1);
