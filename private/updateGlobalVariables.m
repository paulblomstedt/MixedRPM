function updateGlobalVariables(inds, i2, priorPar, data, cdat, nCategories, nAttributes)
% Updates global variables when items inds are moved into cluster i2.
% Before moving, the items to be moved have to have been in the same cluster.
% inds needs to be a row vector.

% Author(s): Paul Blomstedt, Pekka Marttinen

global PARTITION; 
global COUNTS; 
global SUMCOUNTS;
global SAMPLEVAR; 
global SAMPLEMEAN; 
global ADDITION_DIFFERENCE;
global REMOVAL_DIFFERENCE;
global JOIN_DIFFERENCE;

i1 = PARTITION(inds(1));
PARTITION(inds)=i2;
inds_i1 = find(PARTITION==i1); 
inds_i2 = find(PARTITION==i2); 

diffInCounts = computeDiffInCounts(inds, size(COUNTS,2), data, nCategories);
diffInSumCounts = sum(diffInCounts);

COUNTS(:,:,i1) = COUNTS(:,:,i1) - diffInCounts;
SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:) - diffInSumCounts;
COUNTS(:,:,i2) = COUNTS(:,:,i2) + diffInCounts;
SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:) + diffInSumCounts;
[SAMPLEMEAN(i1,:), SAMPLEVAR(i1,:)] = calcStats(inds_i1,cdat); 
[SAMPLEMEAN(i2,:), SAMPLEVAR(i2,:)] = calcStats(inds_i2,cdat); 

updateLogmlTable([i1 i2], priorPar, nCategories, nAttributes);

REMOVAL_DIFFERENCE(PARTITION==i1) = nan;
REMOVAL_DIFFERENCE(PARTITION==i2) = nan;
ADDITION_DIFFERENCE(:,[i1 i2]) = nan;

npops = size(COUNTS,3); 
partind = unique(PARTITION); 
nonempty = find(ismember(1:npops,partind)); %nonemtpy clusters
JOIN_DIFFERENCE(nonempty,i2) = nan; 
JOIN_DIFFERENCE(i2,nonempty) = nan; 
if isempty(find(PARTITION==i1, 1))
    % i1 became empty
    JOIN_DIFFERENCE(:,i1) = 0;
    JOIN_DIFFERENCE(i1,:) = 0;
    JOIN_DIFFERENCE(i1,i1) = nan;
else
    JOIN_DIFFERENCE(:,i1) = nan;
    JOIN_DIFFERENCE(i1,:) = nan;
end