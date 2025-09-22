function [OutPut,r,p]=LuPairRegressPlot3D_Group(data1,data2,dataGroup,varargin)
% LuPairRegressPlot3D_Group: Multiple regression with 3D scatter and plane fit
%
%   data1: [n x 2] predictor variables
%   data2: [n x 1] response variable
%   dataGroup: grouping vector for scatter colors
%   varargin: {Param}, {Param,Cov}

% Ensure correct shape
data1 = data1(:,:);
data2 = data2(:);

% -------------------- Default Parameters --------------------
if nargin==3
   Param.Color=[0.1 0.1 0.1];
   Param.Marker='o';
   Param.MarkerSize=20;
   Param.Rtype='pearson';
   Param.xLabel='X1';
   Param.yLabel='X2';
   Param.zLabel='Y';
   Param.View=[45 30]; % default 3D view
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

% -------------------- Scatter Plot by Group --------------------
 hold on;
uniqueGroups = unique(dataGroup);
for iGroup=1:length(uniqueGroups)
    temp1=find(dataGroup==uniqueGroups(iGroup));
    if ~isempty(temp1)
        scatter3(data1(temp1,1),data1(temp1,2),data2(temp1),Param.MarkerSize,...
            'MarkerEdgeColor',Param.Color(iGroup,:),...
            'MarkerFaceColor',Param.Color(iGroup,:),...
            'Marker',Param.Marker);
    end
end
xlabel(Param.xLabel);
ylabel(Param.yLabel);
zlabel(Param.zLabel);
grid on;

% -------------------- Multiple Linear Regression --------------------
[B,BINT,R,RINT,STATS]=regress(data2,[ones(size(data1,1),1) data1]);
OutPut.B=B;
OutPut.BINT=BINT;
OutPut.R=R;
OutPut.RINT=RINT;
OutPut.STATS=STATS;

% -------------------- Regression Plane --------------------
xRange = linspace(min(data1(:,1)),max(data1(:,1)),20);
yRange = linspace(min(data1(:,2)),max(data1(:,2)),20);
[X,Y] = meshgrid(xRange,yRange);
Z = B(1) + B(2)*X + B(3)*Y;
mesh(X,Y,Z,'FaceAlpha',0.2,'EdgeColor','none','FaceColor',[0.5 0.5 0.5]);

% -------------------- Correlation --------------------
predY = [ones(size(data1,1),1) data1]*B;

[r,p]=corr(data2,predY,'rows','complete','type',Param.Rtype);

if ~isempty(Cov)
    [r(end+1),p(end+1)] = partialcorr(data2,data1(:,1),[data1(:,2) Cov],'rows','complete','type',Param.Rtype);
    [r(end+1),p(end+1)] = partialcorr(data2,data1(:,2),[data1(:,1) Cov],'rows','complete','type',Param.Rtype);
else
    [r(end+1),p(end+1)] = partialcorr(data2,data1(:,1),[data1(:,2)],'rows','complete','type',Param.Rtype);
    [r(end+1),p(end+1)] = partialcorr(data2,data1(:,2),[data1(:,1)],'rows','complete','type',Param.Rtype);

end

% -------------------- Annotation --------------------

title(['Overall fit: r = ' showNum(r(1),3) ', p' showPvalue(p(1),3)],'Color','k');
grid on;
view(Param.View);

if length(p)>1

    strX1 =['Par. r(Z'  ',' Param.xLabel ') = ' showNum(r(2),3) ', p' showPvalue(p(2),3)];
    strX2 = ['Par. r(Z'  ',' Param.yLabel ') = ' showNum(r(3),3) ', p' showPvalue(p(3),3)];

   if p(2)<0.05 
      if r(2)>0
      text(min(xRange),min(yRange),max(data2),...
        strX1,'FontSize',10,'Color','r','HorizontalAlignment','left','HorizontalAlignment','left');
      else
       text(min(xRange),min(yRange),max(data2),...
        strX1,'FontSize',10,'Color','b','HorizontalAlignment','left','HorizontalAlignment','left');
         
      end
   else
        text(max(xRange),min(yRange),max(data2),strX1,'FontSize',10,'Color','k','HorizontalAlignment','left');

   end
   if p(3)<0.05 
      if r(3)>0
      text(min(xRange),max(yRange),max(data2),strX2,'FontSize',10,'Color','r','HorizontalAlignment','left','VerticalAlignment','top');
      else
      text(min(xRange),max(yRange),max(data2),strX2,'FontSize',10,'Color','b','HorizontalAlignment','left');
         
      end
   else
      text(min(xRange),max(yRange),max(data2),strX2,'FontSize',10,'Color','k','HorizontalAlignment','left');


   end
end

end

% function txt = oneLineText(str)
% % Utility to replace newlines or unwanted breaks with spaces
%     txt = regexprep(str,'[\n\r]+',' ');
% end
