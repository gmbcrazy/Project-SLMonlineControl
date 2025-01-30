function TimeBarPlot(Data,Color,Xstart,Ystart,Xwidth,Ywidth, varargin)

if nargin==7
   AddArrow=1;
   ArrowVec=varargin{1};
   AddSham=0;
elseif nargin==8
   AddArrow=1;
   ArrowVec=varargin{1};
   AddSham=1;
   ShamVec=varargin{2};


else
   AddArrow=0;
   AddSham=0;

end

if size(Color,1)==1
   Color=repmat(Color,length(Xwidth),1);
   ZeroColor=[0.9 0.9 0.9];
else
   ZeroColor=[0.9 0.9 0.9];
end

StartAll=[];
for i=1:length(Data)
    StartAll=[StartAll;Xstart];
    if Data(i)~=0
       barplot(Xstart,Ystart,Xwidth(i),Ywidth*0.6,ZeroColor,1);hold on;

       if AddSham==0
          plot([Xstart Xstart],[Ystart Ystart+Ywidth*0.6],'Color',Color(Data(i),:),'LineWidth',2);
       else
          if ShamVec(i)==0
             plot([Xstart Xstart],[Ystart Ystart+Ywidth*0.6],'Color',Color(Data(i),:),'LineWidth',2);
          else
             plot([Xstart Xstart],[Ystart Ystart+Ywidth*0.6],'Color',[Color(Data(i),:) 0.2],'LineWidth',2);
         
          end
       end

    else
       barplot(Xstart,Ystart,Xwidth(i),Ywidth*0.6,ZeroColor,1);hold on;    
    end
    Xstart=Xstart+Xwidth(i);
end





Xtemp=cumsum(Xwidth(:));
set(gca,'xtick',[0;Xtemp],'tickdir','out','ytick',[],'xlim',[0 Xtemp(end)],'ylim',[Ystart, Ywidth+Ystart],'box','off');

ax = gca;
if AddArrow==1
   I1=find(ArrowVec>0);
   for i=1:length(I1)

        xPos = StartAll(I1(i)); % Position at the end of the bar
        yPos = Ywidth+Ystart; % Sequence position in the plot
        xNorm = ax.Position(1) + ax.Position(3) * (xPos - ax.XLim(1)) / diff(ax.XLim);
        yNorm = ax.Position(2) + ax.Position(4) * (yPos - ax.YLim(1)) / diff(ax.YLim);
        
        % Create the annotation
        % annotation('textarrow', [xNorm, xNorm], [yNorm, yNorm - 0.05], 'String', 'VolOut', 'FontSize', 8);
        annotation('textarrow', [xNorm, xNorm], [yNorm, yNorm - 0.3], 'FontSize', 6);

   end


end




ax.YAxis.Visible = 'off';
end



function barplot(Xstart,Ystart,Xwidth,Ywidth,color,alpha)


for i=1:length(Xstart)
    fill([Xstart(i) Xstart(i)+Xwidth(i) Xstart(i)+Xwidth(i) Xstart(i) Xstart(i)],[Ystart(i) Ystart(i) Ystart(i)+Ywidth(i) Ystart(i)+Ywidth(i) Ystart(i)],color);
    h=get(gca,'children');
    set(h(1),'edgecolor',color,'edgealpha',0,'facealpha',alpha);
hold on;
end
end
% fill([X(1) X(1:end) fliplr([X(1:end) X(end)])],[Yabove(1) Yabove fliplr([Ylow Ylow(end)])],color);
% hold on;
% plot(X,Y,'color',color);
% h=get(gca,'children');
% set(h(2),'facealpha',alpha,'linestyle','-');
% set(h(2),'edgecolor',color,'edgealpha',alpha);
% set(h(1),'linestyle','-');
% 



