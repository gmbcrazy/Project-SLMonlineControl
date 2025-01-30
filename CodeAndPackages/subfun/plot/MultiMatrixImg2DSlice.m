function MultiMatrixImg2DSlice(Data,Xstep)
%%%%%%%%Data is 3D matrix, each sample is Data(:,:,i)



%# plot each slice as a texture-mapped surface (stacked along the Z-dimension)
Xstart=0;
hold on;
for k=1:length(Data)
    Xtemp=Xstart+[1:size(Data{k},1)];
    Xstart=Xtemp(end)+Xstep;
    Ytemp=1:size(Data{k},2);
    imagesc(Xtemp,Ytemp,Data{k});
    hold on;
end
