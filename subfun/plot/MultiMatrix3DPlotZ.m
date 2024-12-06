function MultiMatrix3DPlotZ(Data,ZPlot,Alpha)
%%%%%%%%Data is 3D matrix, each sample is Data(:,:,i)

[X,Y]=meshgrid(1:size(Data,1),1:size(Data,2));

% surf(x,y,[1 2;3 4],gradient(z))
%# create stacked images (I am simply repeating the same image 5 times)
% img = load('clown');
I = repmat(squeeze(Data(:,:,1)),[1 1 size(Data,3)]);
% cmap = img.map;

%# coordinates
[X,Y] = meshgrid(1:size(I,2), 1:size(I,1));
Z = zeros(size(I,1),size(I,2));

%# plot each slice as a texture-mapped surface (stacked along the Z-dimension)
for k=1:size(I,3)
    % surface('XData',X, 'YData',Y, 'ZData',Z+ZPlot(k), ...
        % 'CData',Data(:,:,k), 'CDataMapping','direct', ...
        % 'EdgeColor','none', 'FaceColor','texturemap','FaceAlpha',1);
        surface('XData',X, 'YData',Y, 'ZData',Z+ZPlot(k), ...
        'CData',Data(:,:,k), 'CDataMapping','scaled', ...
        'EdgeColor','none','FaceAlpha',Alpha);

    % surface('XData',floor(Z.*XPlot(k)*size(Z,1)), 'YData',X, 'ZData',Y, ...
    %     'CData',Data(:,:,k)*size(Data,1)/2, 'CDataMapping','direct', ...
    %     'EdgeColor','none','FaceColor','texturemap','FaceAlpha',0.3);

    hold on;
end
% colormap('parula')
% view(3), box off,
% axis tight square
% set(gca, 'YDir','reverse', 'ZLim',[1 size(Data,3)])
