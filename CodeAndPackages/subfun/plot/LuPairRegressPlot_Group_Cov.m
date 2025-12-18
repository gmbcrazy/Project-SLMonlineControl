function [OutPut,r,p,h]=LuPairRegressPlot_Group_Cov(data1,data2,Cov,dataGroup,varargin)    

data1=data1(:);
data2=data2(:);

if nargin==4
   Param.Color=[0.1 0.1 0.1];
   Param.Marker='c';
   Param.MarkerSize=10;
   Param.Rtype='pearson';
   Param.xLim=[min(data1) max(data1)];
   Param.yLim=[min(data2) max(data2)];
   Param.xLabel=[];
   Param.yLabel=[];
   % Cov=[];
   Param.RawcorrPlot=1;

elseif nargin==5
   Param=varargin{1};
   % Cov=varargin{2};

else

end

if ~isfield(Param,'RawcorrPlot')
   Param.RawcorrPlot=1;
end
if isempty(Param.xLim)
   AddIn=(max(data1)-min(data1))/10;
    Param.xLim=[min(data1)-AddIn max(data1)+AddIn];
end
if isempty(Param.yLim)
   AddIn=(max(data2)-min(data2))/10;

    Param.yLim=[min(data2)-AddIn max(data2)+AddIn];
end


[B,BINT,R,RINT,STATS]=regress(data2,[ones(length(Cov),1) Cov]);

if Param.RawcorrPlot==1

   Param.xLim=[];
   Param.yLim=[];

h(1)=subplot(1,2,1);
[OutPut,r,p]=LuPairRegressPlot_Group(data1,data2,dataGroup,Param)
h(2)=subplot(1,2,2);
[B,BINT,R,RINT,STATS]=regress(data2,[ones(length(Cov),1) Cov]);

% Param.yLim=[min(R) max(R)];
[OutPut,r,p]=LuPairRegressPlot_Group(data1,R,dataGroup,Param);
else
   Param.xLim=[];
   Param.yLim=[];

[OutPut,r,p]=LuPairRegressPlot_Group(data1,R,dataGroup,Param);

h=NaN;
end
