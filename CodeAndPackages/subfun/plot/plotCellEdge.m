function [Ix,Iy]=plotCellEdge(imC,cellIDmark,varargin)



if nargin<3
   colorCell=jet(length(cellIDmark));
   Markersize=3;
   LabelCell=0;
elseif nargin==3
   colorCell=varargin{1};
   Markersize=3;
   LabelCell=0;
elseif nargin==4
   colorCell=varargin{1};
   Markersize=varargin{2};
   LabelCell=0;
elseif nargin==5
   colorCell=varargin{1};
   Markersize=varargin{2};
   LabelCell=varargin{3};

else

end
hold on;

for i=1:length(cellIDmark)
    Ix{i}=[];
    Iy{i}=[];
    [Ix{i},Iy{i},~]=find(imC==cellIDmark(i));
    if ~isempty(Ix{i})
       % if LabelCell==0
       plot(Iy{i}+1,Ix{i}+1,'.','color',colorCell(i,:),'Markersize',3);
       % end
       % plot(median(Iy{i}+1),median(Ix{i}+1),'o','color',colorCell(i,:),'Markersize',6);

       if LabelCell~=0
          text(median(Iy{i}),median(Ix{i}),num2str(cellIDmark(i)),'FontSize',12,'Color',colorCell(i,:),'FontWeight','bold');
       end
    end
end
hold off;

