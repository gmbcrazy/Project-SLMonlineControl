
function [XYPosPixel,Z]=gplXYtoPixel(tbl,yaml)

% Extract X, Y, and Z coordinates from Pos3D
% Extract relevant parameters from yaml structure
% Convert full field into imaging FOV

ScanAmp_X = yaml.ScanAmp_X;
ScanAmp_Y = yaml.ScanAmp_Y;
ScanAmp_V_FOV_X = (ScanAmp_X - mean(ScanAmp_X)) + mean(ScanAmp_X);  % Centre, scale, offset
ScanAmp_V_FOV_Y = (ScanAmp_Y - mean(ScanAmp_Y)) + mean(ScanAmp_Y);  % Centre, scale, offset

FOVsize_OpticalZoom = yaml.FOVsize_OpticalZoom;
FOVsize_PX = yaml.FOVsize_PX;
FOVsize_UM_1x = yaml.FOVsize_UM;
% Build Look-Up Tables (LUTs)
LUTx = linspace(ScanAmp_V_FOV_X(1), ScanAmp_V_FOV_X(2), FOVsize_PX);
LUTy = linspace(ScanAmp_V_FOV_Y(1), ScanAmp_V_FOV_Y(2), FOVsize_PX);


Xv=tbl.X;
Yv=tbl.Y;

Z =[];
if ismember('Z', tbl.Properties.VariableNames)
   Z = tbl.Z;
   Zplane=unique(tbl.Z);
   [~,Zplane]=ismember(Z,Zplane);
   Z = [Z(:) Zplane];
end


% Convert full field into imaging FOV
ScanAmp_V_FOV_X = (ScanAmp_X - mean(ScanAmp_X)) + mean(ScanAmp_X);  % Centre, scale, offset
ScanAmp_V_FOV_Y = (ScanAmp_Y - mean(ScanAmp_Y)) + mean(ScanAmp_Y);  % Centre, scale, offset



% Convert pixel coordinates to voltages
[~, Xpx] = histc(Xv, LUTx);
[~, Ypx] = histc(Yv, LUTy);

XYPosPixel=[Xpx(:) Ypx(:)];
