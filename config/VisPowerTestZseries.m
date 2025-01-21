
Data=0:10;
Xstart=[]

TimeBarPlot(Data,Color,Xstart,Ystart,Xwidth,Ywidth, varargin)


SLMnum=10;
Xstart=1;
Xwidth=repmat(31,1,SLMnum+1);
Ystart=0;
Ywidth=1;
Color=[0 0 0];

TimeBarPlot([0;repmat(1,SLMnum,1)],Color,1,Ystart,Xwidth,Ywidth)

