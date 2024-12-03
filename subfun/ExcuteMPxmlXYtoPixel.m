
function XYPosPixel=ExcuteMPxmlXYtoPixel(xy,confSet)

% Extract X, Y, and Z coordinates from Pos3D
% Extract relevant parameters from yaml structure
% Convert full field into imaging FOV

if isstruct(confSet)
numX = confSet.SLM_Pixels_X;
numY = confSet.SLM_Pixels_Y;
elseif isnumeric(confSet)
numX = confSet(1);
numY = confSet(2);

else

end
if istable(xy)
Xpx=[ceil(xy.X*numX) ];
Ypx=[ceil(xy.Y*numY) ];
    
elseif isnumeric(xy)
Xpx=[ceil(xy(:,1)*numX) ];
Ypx=[ceil(xy(:,2)*numY) ];

else

end
XYPosPixel=[Xpx(:) Ypx(:)];
