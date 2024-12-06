function [im,imC]=pltImageCell(stat, ops, iscell,cellID)


 validC=find(iscell(:,1)>0);
 nCells=length(validC);

if nargin<4
   % cellID(validC)=1:validC;
   iscell(validC,1)=1:nCells;
end
    meanImgE = ops.meanImgE;
%    imagesc([1:ops.Lx], [1:ops.Ly], meanImgE);
    im = zeros(ops.Ly, ops.Lx);
    % ncells = numel(stat);
    countValid=1;
    for n = 1:nCells
        % if iscell(n, 1) > 0.5
            ypix = stat{validC(n)}.ypix(~stat{validC(n)}.overlap);
            xpix = stat{validC(n)}.xpix(~stat{validC(n)}.overlap);
            % im(ops.Ly - ypix, xpix) = countValid;
            % im(ops.Ly - ypix, xpix) = cellID(n);

            % im(ops.Ly - ypix, xpix) = countValid;
            for ipix=1:length(ypix)
            im(ypix(ipix), xpix(ipix)) = cellID(validC(n));
            end
        % end
    end
    % llx=size(im,1);
    % lly=size(im,2);
    % imTemp=zeros(llx+1,lly+1);
    [dZdx, dZdy] = gradient(im);
    % imTempdx=imTemp;
    % imTempdy=imTemp;
    % % % imTempdx(2:llx+1,2:lly+1)=dZdx;
    % % % imTempdy(2:llx+1,2:lly+1)=dZdy;
    % imTempdx(1:llx,1:lly)=dZdx;
    % imTempdy(1:llx,1:lly)=dZdy;
    % 
    % dZdx=imTempdx(2:llx+1,2:lly+1);
    % dZdy=imTempdy(2:llx+1,2:lly+1);


    imC=im;
    imC(abs(dZdx)==0&abs(dZdy)==0)=0;
end



