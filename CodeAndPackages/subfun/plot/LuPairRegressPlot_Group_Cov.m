function [OutPut,r,p,h]=LuPairRegressPlot_Group_Cov(data1,data2,Cov,dataGroup,varargin)    

data1=data1(:);
data2=data2(:);

if nargin==4
   Param.Color=[0.1 0.1 0.1];
   Param.Marker='s';
   Param.MarkerSize=2;
   Param.Rtype='pearson';
   Param.xLim=[min(data1) max(data1)];
   Param.yLim=[min(data2) max(data2)];
   Param.xLabel=[];
   Param.yLabel=[];
   % Cov=[];


elseif nargin==5
   Param=varargin{1};
   % Cov=varargin{2};

else

end

h(1)=subplot(1,2,1);
[OutPut,r,p]=LuPairRegressPlot_Group(data1,data2,dataGroup,Param)
h(2)=subplot(1,2,2);
[B,BINT,R,RINT,STATS]=regress(data2,[ones(length(Cov),1) Cov]);
[OutPut,r,p]=LuPairRegressPlot_Group(data1,R,dataGroup,Param);


