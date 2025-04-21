function OutTBL=MPSeqFolder_1TargetXNon(Folder,varargin)

%%Get MarkPoints Seq Info from a given Tseries folder

if nargin==1
   xyPixelNum=[512;512];
elseif nargin==2
   xyPixelNum=varargin{1};  %xyPixelNum could be a configuration variable read from ReadYaml('SLMsetting.yml');
   if isempty(xyPixelNum)
      xyPixelNum=[512;512];
   end
elseif nargin==3
   xyPixelNum=varargin{1};
   AllTestPoints3D=varargin{2};
else

end

%%Get SLM excute information, MarkPoints group include 1 target and X nontargets.

MPxmlFiles=dir([Folder '*MarkPoints*.xml']);

for i = 1:length(MPxmlFiles)
        % Extract the cycle number from the MarkPoints.xml file name
        xmlFileName = MPxmlFiles(i).name;
        xmlTokens = regexp(xmlFileName, 'Cycle(\d{5})', 'tokens');
        if ~isempty(xmlTokens)
            markCycle = str2double(xmlTokens{1}{1});
            [tbl,StimID{i},StimPowerTemp]=MPxml2Table([MPxmlFiles(i).folder '\' MPxmlFiles(i).name]);
            XYPosPixel=ExcuteMPxmlXYtoPixel([tbl.X(1) tbl.Y(1)], xyPixelNum);        %The 1st MarkPoints is the cell targets.
            ExcuteMP(i,:)=[markCycle XYPosPixel str2num(StimPowerTemp)]; 
        end
        % StimuliPower=[StimuliPower;MPxml2yaml([MPxmlFiles(i).folder '\' MPxmlFiles(i).name])];
end

if exist('AllTestPoints3D')
   for i = 1:length(MPxmlFiles)
        temp1=sum(abs(ExcuteMP(i,2:3)-AllTestPoints3D(:,1:2)).^0.5,2);

        [temp2,I1]=min(temp1);
%         if temp2<=4                         %%<4 pixels difference of xy coordinates compared to original MarkPoints
        if temp2<=5                         %%<4 pixels difference of xy coordinates compared to original MarkPoints

           ExcuteMP(i,5)=I1;
        else
           ExcuteMP(i,5)=0;
        end
   end
else
   ExcuteMP(i,5)=0;
end

OutTBL=array2table(ExcuteMP,'VariableNames',{'markCycle','X','Y','UncagingLaserPower','Point'});
