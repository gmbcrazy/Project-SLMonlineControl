function [map_idx,overlap] = mapCellStat1_To_CellStat2(CellStat1, CellStat2)
% map_idx(j) = index i in CellStat1 that matches CellStat2{j}
% map_idx(j) = 0 if CellStat2{j} is empty or no match

n2 = numel(CellStat2);
n1 = numel(CellStat1);

map_idx = zeros(n2,1);   % default = 0
overlap = map_idx;

for j = 1:n2
    c = CellStat2{j};

    % empty entry -> leave as 0
    if isempty(c)
        continue
    end

    for i = 1:n1
        s = CellStat1{i};
        if isempty(s)
            continue
        end


        % overlap of x pixels
        overlap(j) = min(numel(intersect(c.xpix, s.xpix)) / numel(union(c.xpix, s.xpix)),...
            numel(intersect(c.ypix, s.ypix)) / numel(union(c.ypix, s.ypix)));



        if overlap(j) > 0.99   % threshold; adjust as needed
            map_idx(j) = i;
            break          % found match → go to next j
        end
    end
end
end
