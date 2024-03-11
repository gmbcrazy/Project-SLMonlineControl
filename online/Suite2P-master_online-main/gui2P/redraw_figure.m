
function I = redraw_figure(h)
% 
Sat1     =  ones(h.dat.cl.Ly, h.dat.cl.Lx);
Sat2     =  ones(h.dat.cl.Ly, h.dat.cl.Lx);
H1              = zeros(h.dat.cl.Ly, h.dat.cl.Lx);
H2              = zeros(h.dat.cl.Ly, h.dat.cl.Lx);

[iclust1, iclust2, V1, V2] = ...
    getviclust(h.dat.stat, h.dat.cl.Ly,  h.dat.cl.Lx, h.dat.cl.vmap, h.dat.F.ichosen);
%%%

%%%%
iselect     = iclust1==h.dat.F.ichosen;
Sat1(iselect)= 0;

iselect     = iclust2==h.dat.F.ichosen;
Sat2(iselect)= 0;


H1(iclust1>0)   = h.dat.cl.rands(iclust1(iclust1>0));
H2(iclust2>0)   = h.dat.cl.rands(iclust2(iclust2>0));
%%%
 H1=flip(H1,2);
Sat1=flip(Sat1,2);
V1=flip(V1,2);
 H2=flip(H2,2);
Sat2=flip(Sat2,2);
V2=flip(V2,2);
%%%%
I = hsv2rgb(cat(3, H1, Sat1, V1));
I = min(I, 1);
axes(h.axes2); imagesc(I);
xlim([h.dat.xlim]); ylim([h.dat.ylim]);
axis off
%%%%%modified by Hari 06.10.2022, comment by Lu
% figure;
% imshow(mean(I,3),[])

Plane_fileIscell=extractfield(h.dat.stat,'iscell');
for is=1:size(Plane_fileIscell,2)
    Plane_fileIscellTemp=cell2mat(Plane_fileIscell(is));
    Plane_fileIscellTotal(is)=Plane_fileIscellTemp;
end
Plane_fileIscell=Plane_fileIscellTotal;

Plane_fileXYcorr=extractfield(h.dat.stat,'med');
Plane_fileXYcorrX=Plane_fileXYcorr(1:2:end);
Plane_fileXYcorrY=Plane_fileXYcorr(2:2:end);
Plane_fileXYcorr=[Plane_fileXYcorrX(:)  Plane_fileXYcorrY(:)];
%Plane_filePb=extractfield(h.dat.stat,'mrs'); 
  icount=1;
for i=1:size(Plane_fileIscell,2)
Plane_fileIscelltemp=Plane_fileIscell(i);
%Plane_filePbtemp=Plane_filePb(i);
     if Plane_fileIscelltemp>0
       Plane_fileXYcorrtemp=Plane_fileXYcorr(i,:);
       Plane_fileXYcorrTotal(icount,:)=Plane_fileXYcorrtemp;

       %Plane_fileXYcorrTotal(:,2)=[512-Plane_fileXYcorrTotal(:,2)];
       icount=icount+1;
     %  imshow(Plane_file.ops.sdmov,[])
       hold on
      %plot(Plane_fileXYcorrtemp(2),Plane_fileXYcorrtemp(1),'Or')
      plot(512-Plane_fileXYcorrtemp(2),Plane_fileXYcorrtemp(1),'Or')%I need to flip it 
    end
end
%% Comment by Lu
% figure;
% imshow(V1,[])
% roi_img=imshow(V1,[]);
% imwrite(roi_img.CData,[h.dat.ops.ResultsSavePath,'mask_img_Plane.tif'])
%% Comment by Lu

%%for Naparm points
Plane_fileXYcorrTotal=flip(Plane_fileXYcorrTotal,2);
Plane_fileXYcorrTotal(:,1)=512-Plane_fileXYcorrTotal(:,1);
Plane_fileXYcorrTotal=round(Plane_fileXYcorrTotal);
points.h=[];
points.X=(Plane_fileXYcorrTotal(:,1))';
points.Y=(Plane_fileXYcorrTotal(:,2))';
points.Z=ones(1,size(Plane_fileXYcorrTotal,1));
points.Zum=99.9*(ones(1,size(Plane_fileXYcorrTotal,1)));
points.OffsetX=ones(1,size(Plane_fileXYcorrTotal,1));
points.OffsetY=ones(1,size(Plane_fileXYcorrTotal,1));
points.Idx=1:1:size(Plane_fileXYcorrTotal,1);
points.Img=ones(1,size(Plane_fileXYcorrTotal,1));
points.Group=nan(1,size(Plane_fileXYcorrTotal,1));
points.GroupCentroidX=nan(1,size(Plane_fileXYcorrTotal,1));
points.GroupCentroidY=nan(1,size(Plane_fileXYcorrTotal,1));
points.Counter=size(Plane_fileXYcorrTotal,1);
points.Weight=ones(1,size(Plane_fileXYcorrTotal,1));
points.Selected=zeros(1,size(Plane_fileXYcorrTotal,1));
%%Here Hari generate random point in each group for testing purpose 
roi_size=numel(points.Group);
roi_idx=1:roi_size;
if roi_size<50
roi_sort=1:roi_size;
roi_rand=roi_sort(sort(randperm(roi_size,5)));
roi_rand_diff=diff(roi_rand);
roi_rand_diff(end+1)=roi_size-sum(roi_rand_diff);
sum(roi_rand_diff);

elseif  roi_size>50 && roi_size<100
    roi_sort=1:roi_size;
roi_rand=roi_sort(sort(randperm(roi_size,12)));
roi_rand_diff=diff(roi_rand);
roi_rand_diff(end+1)=roi_size-sum(roi_rand_diff);
sum(roi_rand_diff);

else 
    roi_sort=1:roi_size;
    roi_rand=roi_sort(sort(randperm(roi_size,20)));
    roi_rand_diff=diff(roi_rand);
    roi_rand_diff(end+1)=roi_size-sum(roi_rand_diff);
    sum(roi_rand_diff);
end

roi_rand_diffCumsum=(cumsum(roi_rand_diff))+1;
for g=1:numel(roi_rand_diff)
    roi_sizeTemp=roi_idx(roi_rand_diffCumsum(g)-roi_rand_diff(g):roi_rand_diffCumsum(g)-1); 
    roi_group(roi_sizeTemp)=g;
end
points.Group=roi_group;%% to generate the groups of different points

save([h.dat.ops.ResultsSavePath,'POINT.mat'],'points')
points={'X','Y','Z','Zum','OffsetX','OffsetY','Idx','Img','Group','GroupCentroidX','GroupCentroidY','Counter','Weight','Selected'};
%%%%%%%%%%%%%%%%%%%%%%%%%%Hari, comment by Lu
I = hsv2rgb(cat(3, H2, Sat2, V2));
I = min(I, 1);
axes(h.axes3); imagesc(I);
xlim([h.dat.xlim]); ylim([h.dat.ylim]);
axis off

end

function [iclust1, iclust2, V1, V2] = getviclust(stat, Ly, Lx, vmap, ichosen)

iclust1 = zeros(Ly, Lx);
iclust2 = zeros(Ly, Lx);
V1      = zeros(Ly, Lx);
V2      = zeros(Ly, Lx);

for j = 1:numel(stat)
    ipix    = stat(j).ipix;
    lambda   = stat(j).lambda;
    
    if ichosen==j
        inew = true(numel(ipix), 1);
    else
        if stat(j).iscell
            inew    = lambda(:)>(V1(ipix) + 1e-6);
        else
            inew    = lambda(:)>(V2(ipix) + 1e-6);
        end
    end
    
    switch vmap
        case 'var'
            L0      = stat(j).lambda(inew);
        case 'unit'
            L0      = stat(j).lam(inew);    
    end
    if stat(j).iscell
        V1(ipix(inew))      = L0;
        iclust1(ipix(inew)) = j;
    else
        V2(ipix(inew))      = L0;
        iclust2(ipix(inew)) = j;
    end
end
mV = mean([V1(V1>0); V2(V2>0)]);
V1 = V1/mV;
V2 = V2/mV;
% iclust1=flip(iclust1,2);
% iclust2=flip(iclust2,2);
% V1=flip(V1,2);
% V2=flip(V2,2);
end
%%

