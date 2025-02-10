function MultibarPlot(Data,Color,Xstart,Ystart,Xwidth,Ywidth)

for i=1:length(Data)
    if Data(i)~=0
       barplot(Xstart+(i-1)*Xwidth,Ystart,Xwidth,Ywidth,Color(Data(i),:),1);hold on;
    else
       barplot(Xstart+(i-1)*Xwidth,Ystart,Xwidth,Ywidth,[0.7 0.7 0.7],1);hold on;    
    end
end