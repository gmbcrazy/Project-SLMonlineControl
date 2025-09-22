function [OutPut,r,p]=LuPairRegressPlot_Group(data1,data2,dataGroup,varargin)    

% Ensure column vectors
data1=data1(:);
data2=data2(:);

% -------------------- Default Parameters --------------------
if nargin==3
   Param.Color=[0.1 0.1 0.1];
   Param.Marker='s';
   Param.MarkerSize=2;
   Param.Rtype='pearson';
   Param.xLim=[min(data1) max(data1)];
   Param.yLim=[min(data2) max(data2)];
   Param.xLabel=[];
   Param.yLabel=[];
   Cov=[];
elseif nargin==4
   Param=varargin{1};
   Cov=[];
elseif nargin==5
   Param=varargin{1};
   Cov=varargin{2};
else
   error('Invalid number of input arguments');
end

% -------------------- Axis Limits --------------------
if isempty(Param.xLim)
    Param.xLim=[min(data1) max(data1)];
end
if isempty(Param.yLim)
    Param.yLim=[min(data2) max(data2)];
end

% -------------------- Scatter Plot by Group --------------------
uniqueGroups = unique(dataGroup);
for iGroup=1:length(uniqueGroups)
    temp1=find(dataGroup==uniqueGroups(iGroup));
    if ~isempty(temp1)
       h=scatter(data1(temp1),data2(temp1),Param.MarkerSize,...
           'MarkerEdgeColor',Param.Color(iGroup,:),...
           'MarkerFaceColor',Param.Color(iGroup,:));
       set(h, 'Marker', Param.Marker); 
    end
    hold on;
end

% -------------------- Regression --------------------
[B,BINT,R,RINT,STATS]=regress(data2,[ones(length(data1),1) data1]);
OutPut.B=B;
OutPut.BINT=BINT;
OutPut.R=R;
OutPut.RINT=RINT;
OutPut.STATS=STATS;

x=[min(data1(:)) max(data1(:))]';
y=[[1 1]' x]*B;

% -------------------- Correlation --------------------
[r,p]=corr(data1,data2,'rows','complete','type',Param.Rtype);
if ~isempty(Cov)
   [r(end+1),p(end+1)]=partialcorr(data1,data2,Cov,'rows','complete','type',Param.Rtype);
end

% -------------------- Plotting Correlation --------------------
if r(1)>0
   PColor=[1 0 0]; % Red for positive
elseif r(1)<0
   PColor=[0 0 1]; % Blue for negative
else
   PColor=[0.01 0.01 0.01];
end

if p(1) < 0.05
    text(Param.xLim(2),Param.yLim(2),['r = ' showNum(r(1),3) ', p' showPvalue(p(1),3)],...
        'color',PColor,'fontsize',10,'horizontalalignment','right','verticalalignment','bottom','fontname','Times New Roman');
    hold on;
    plot(x,y,'-','color',PColor,'linewidth',2); 
else 
    text(Param.xLim(2),Param.yLim(2),['r = ' showNum(r(1),3) ', p' showPvalue(p(1),3)],...
        'color',[0.01 0.01 0.01],'fontsize',10,'horizontalalignment','right','verticalalignment','bottom','fontname','Times New Roman');
end

% -------------------- Partial Correlation --------------------
if length(p)>1
    if p(2) < 0.05
        text(Param.xLim(2),Param.yLim(2),['Par. r = ' showNum(r(2),3) ', p' showPvalue(p(2),3)],...
            'color',PColor,'fontsize',10,'horizontalalignment','right','verticalalignment','top','fontname','Times New Roman');
        hold on;
        plot(x,y,'-','color',PColor,'linewidth',2); 
    else 
        text(Param.xLim(2),Param.yLim(2),['Par. r = ' showNum(r(2),3) ', p' showPvalue(p(2),3)],...
            'color',[0.01 0.01 0.01],'fontsize',10,'horizontalalignment','right','verticalalignment','top','fontname','Times New Roman');
    end
end

% -------------------- Final Axis Settings --------------------
set(gca,'xscale','linear','yscale','linear','xlim',Param.xLim,'ylim',Param.yLim)
xlabel(Param.xLabel);
ylabel(Param.yLabel);
