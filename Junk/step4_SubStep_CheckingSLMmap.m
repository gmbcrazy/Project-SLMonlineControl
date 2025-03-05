ShowCol=6;
ShowRow=round(length(PointAll)/ShowCol)+1;
PlaneZ=confSet.ETL+confSet.scan_Z;

figure;
for iPoint=1:length(PointAll)
    Point=PointAll(iPoint);
    subplot(ShowRow,ShowCol,iPoint)
    I1=find((SLMTrialInfo(:,2)==Point&ismember(SLMTrialInfo(:,3),LaserPower))==1);
    I2=find(ismember(SLMTrialInfo(:,4),RoundCheck)==1);
    I1=intersect(I1,I2);

    if ~isempty(I1)
         if length(I1)==1
         PSTHtemp=squeeze(SLMTrialMap(:,:,:,I1));     
         else
         PSTHtemp=squeeze(mean(SLMTrialMap(:,:,:,I1),4));     
         end
         PlaneI=find(abs(Pos3DNeed(Point,3)-PlaneZ)<1)
         MultPlaneIs2DShow1Plane(PSTHtemp, [], Pos3DNeed(Point,:), [], PlaneZ, PlaneI, [0 1 0], PSTHparam.Clim)
         text(confSet.SLM_Pixels_X/2,0, ['P' num2str(Point) ', n = ' num2str(length(I1))],'color',[0 1 0]);
    end
             set(gca,'xtick',[],'ytick',[])

    if iPoint==1
       text(confSet.SLM_Pixels_X/2,confSet.SLM_Pixels_Y/2,PlaneZ(1)-5,['Point ' num2str(Point) 'Laser' num2str(LaserPower)])
    end
end

colormap(PSTHparam.ColorMap);
papersizePX=[0 0 ShowCol*8 ShowRow*8]
