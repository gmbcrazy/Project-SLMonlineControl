function [NonTargets,NonTargetsPlane]=NonTargetGeneration(SavePath,Pos3DRaw,PosPlane,yaml,confSet,varargin)


% NonAvoidRadius=confSet.SpiralSizeUM*2;
NumNonTargets=confSet.NumNonTarget;

RadiusByPixel=ceil(confSet.SpiralSizeUM/yaml.umPerlPixelX/2);

if nargin==7 %% avoid vessel region for nontarget generating
   MeanImg=varargin{1};
   vesselTh=varargin{2};
   NonTargets = generateNewMarkPoints([Pos3DRaw(:,1:2) PosPlane(:)], confSet.RadiusAvoidParam*RadiusByPixel, length(confSet.ETL), NumNonTargets, RadiusByPixel, [confSet.SLM_Pixels_X confSet.SLM_Pixels_X],MeanImg,vesselTh);
else
   NonTargets = generateNewMarkPoints([Pos3DRaw(:,1:2) PosPlane(:)], confSet.RadiusAvoidParam*RadiusByPixel, length(confSet.ETL), NumNonTargets, RadiusByPixel, [confSet.SLM_Pixels_X confSet.SLM_Pixels_X]);
end


NonTargetsPlane=NonTargets(:,3);
for iplane=1:length(confSet.ETL)
    NonTargets(NonTargetsPlane==iplane,3)=confSet.ETL(iplane)+confSet.scan_Z;
end


[NonTargetsPlane,I1]=sort(NonTargetsPlane);
NonTargets=NonTargets(I1,:);
% NonTargetPath=[SavePath  'NonTargets\'];
% mkdir(NonTargetPath)
% SaveNonTargets=[SavePath 'NonTargets\All']
MarkPoints3D_GPLmaker(NonTargets, yaml, confSet,SavePath,[], 'NonTarget');

save([SavePath '.mat'],'NonTargets','NonTargetsPlane');
