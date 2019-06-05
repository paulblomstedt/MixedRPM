function [dist,Z] = calcEuclidDist(data)

% Author: Pekka Marttinen

ninds = size(data,1);
dist = zeros(nchoosek(ninds,2),1); 
pointer = 1;
for i = 1:ninds-1
    row1 = data(i,:);
    for j = i+1:ninds
        row2 = data(j,:);
        dist(pointer) = sqrt(sum((row1-row2).^2)); 
        pointer = pointer+1;
    end
end
Z = linkage(dist'); 