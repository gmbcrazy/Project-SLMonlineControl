function [Mapping,ValidROI,InvalidROI]=ROIMatchApplyThreshold(Mapping,ROIInfo,CorrTh)

I1=find(ROIInfo.corrMatch>=CorrTh);
I2=find(ROIInfo.corrMatch<CorrTh);

ValidROI=AsignROI(ROIInfo,I1);
InvalidROI=AsignROI(ROIInfo,I2);
Mapping=Mapping(I1,:);

end



function AsignedROI=AsignROI(ROI,Ind)

AsignedROI=ROI;
AsignedROI.FixRoiN=ROI.FixRoiN(:,:,Ind);
AsignedROI.FixRoiB=ROI.FixRoiB(Ind);

AsignedROI.MovingRoiN=ROI.MovingRoiN(:,:,Ind);
AsignedROI.MovingRoiB=ROI.MovingRoiB(Ind);

AsignedROI.MovingRoiRegN=ROI.MovingRoiRegN(:,:,Ind);
AsignedROI.MovingRoiRegB=ROI.MovingRoiRegB(Ind);

AsignedROI.corrMatch=ROI.corrMatch(Ind);
end

