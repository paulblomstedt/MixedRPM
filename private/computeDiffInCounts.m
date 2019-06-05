function diffInCounts = computeDiffInCounts(rowIndexes, nloci, data, nValues)

% Author(s): Eero Linna, Paul Blomstedt

diffInCounts = zeros(nValues, nloci); 

for k=1:nValues
    diffInCounts(k, :) = sum(data(rowIndexes,:)==k-1, 1);
end

