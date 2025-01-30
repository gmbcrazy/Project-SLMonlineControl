function [SLMtarget,SLMtargetCellDist]=SLMtargetMatchCell(SLM3D,NeuroPos3D,DistTh)

for iPoint = 1:size(SLM3D)
    planeI=find(NeuroPos3D(:,3)==SLM3D(iPoint,3));
    Cellxy=NeuroPos3D(planeI,1:2);
    MPxy=SLM3D(iPoint,1:2);
    Dist=sum((abs(Cellxy-MPxy)).^2,2).^0.5;
    [tempDist,Icell]=min(Dist);
    SLMtarget(iPoint)=planeI(Icell);
    SLMtargetCellDist(iPoint)=tempDist;
    clear tempDist Icell
end

SLMtarget(SLMtargetCellDist>DistTh)=-1;     %%Too far away from cell, this target get no valid cells
