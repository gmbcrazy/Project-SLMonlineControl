function MultiMatrixRGB3DPlot(GreenData, RedData, ZPlot, Alpha)
%%%%%%%%%% GreenData & RedData are 3D matrices: (:,:,i) is one plane
%%%%%%%%%% ZPlot defines where to stack the slices
%%%%%%%%%% Alpha sets transparency

[X, Y] = meshgrid(1:size(GreenData, 2), 1:size(GreenData, 1)); % x: columns, y: rows

hold on;

for k = 1:size(GreenData, 3)
    % Extract slice from each channel
    G = mat2gray(GreenData(:, :, k));
    R = mat2gray(RedData(:, :, k));
    B = zeros(size(G)); % No blue channel

    % Construct RGB image
    RGB = cat(3, R, G, B);

    % Plot as texture-mapped surface
    surface('XData', X, 'YData', Y, 'ZData', ZPlot(k)*ones(size(X)), ...
        'CData', RGB, 'CDataMapping', 'direct', ...
        'FaceColor', 'texturemap', 'EdgeColor', 'none', ...
        'FaceAlpha', Alpha);
end

axis tight;
axis square;
view(3);
box on;
xlabel('X');
ylabel('Y');
zlabel('Z');
end
