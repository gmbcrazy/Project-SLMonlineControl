function suite2p_online_proc()
Proc_file=load('C:\Data\Hari\20220706\HS8999\20220706\1\F_HS8999_20220706_plane1_proc.mat');

Plane_file=load('C:\Data\Hari\20220706\HS8999\20220706\1\F_HS8999_20220706_plane1.mat');
%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%
%Proc_file_meanImage
figure;
imshow(Plane_file.ops.sdmov,[])
%Proc_file_meanImage correation image
figure;
imshow(Plane_file.ops.Vcorr,[])
%%
Plane_fileIscell=cell2mat(extractfield(Plane_file.stat,'iscell'));

Plane_fileXYcorr=extractfield(Plane_file.stat,'med');
Plane_fileXYcorrX=Plane_fileXYcorr(1:2:end);
Plane_fileXYcorrY=Plane_fileXYcorr(2:2:end);
Plane_fileXYcorr=[Plane_fileXYcorrX(:)  Plane_fileXYcorrY(:)];
Plane_filePb=extractfield(Plane_file.stat,'mrs');
%close all
figure;
imshow(Plane_file.ops.sdmov,[])
icount=1;
for i=1:size(Plane_file.stat,2)
Plane_fileIscelltemp=Plane_fileIscell(i);
Plane_filePbtemp=Plane_filePb(i);
     if any(Plane_fileIscelltemp) && Plane_filePbtemp>0.5
       Plane_fileXYcorrtemp=Plane_fileXYcorr(icount,:);
       icount=icount+1;
     %  imshow(Plane_file.ops.sdmov,[])
       hold on
      plot(Plane_fileXYcorrtemp(1),Plane_fileXYcorrtemp(2),'Or')
    end
 end
Plane_fileIscell=[Plane_file(1:end).stat(1:end).iscell(1:end)];
x

