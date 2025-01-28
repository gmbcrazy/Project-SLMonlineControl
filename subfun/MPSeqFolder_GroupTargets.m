function [OutTBL,XYTrials,PowerWeight]=MPSeqFolder_GroupTargets(Folder,varargin)

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
   AllFunGroupPoints3D=varargin{2};
else

end

%%Get SLM excute information, MarkPoints group include 1 target and X nontargets.

MPxmlFiles=dir([Folder '*MarkPoints*.xml']);
PowerWeight=[];

for i = 1:length(MPxmlFiles)
        % Extract the cycle number from the MarkPoints.xml file name
        xmlFileName = MPxmlFiles(i).name;
        xmlTokens = regexp(xmlFileName, 'Cycle(\d{5})', 'tokens');
        if ~isempty(xmlTokens)
            markCycle = str2double(xmlTokens{1}{1});
            [tbl,StimID{i},StimPowerTemp,PowerWeightTemp]=MPxml2Table([MPxmlFiles(i).folder '\' MPxmlFiles(i).name]);
            XYPosPixel=ExcuteMPxmlXYtoPixel([tbl.X tbl.Y], xyPixelNum);        %The 1st MarkPoints is the cell targets.
            ExcuteMP(i,:)=[markCycle str2num(StimPowerTemp) NaN]; 
            XYTrials{i}=XYPosPixel;
            PowerWeight(i,:)=str2num(PowerWeightTemp);
        end
        % StimuliPower=[StimuliPower;MPxml2yaml([MPxmlFiles(i).folder '\' MPxmlFiles(i).name])];
end

nPoints=size(AllFunGroupPoints3D{1},2);
FunGroupID=nan(length(MPxmlFiles),1);
if exist('AllFunGroupPoints3D')
   for i = 1:length(MPxmlFiles)
       nPoints=size(XYTrials{i},1);
        for j = 1:length(AllFunGroupPoints3D)
            temp1=sum(abs(XYTrials{i}-AllFunGroupPoints3D{j}(:,1:2)),2);
            if sum(temp1<=4)>nPoints-1
               FunGroupID(i)=j;
               ExcuteMP(i,end)=j;
               break                
            end
        end
   end

end

OutTBL=array2table(ExcuteMP,'VariableNames',{'markCycle','UncagingLaserPower','Group'});
