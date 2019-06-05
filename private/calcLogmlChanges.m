function muutokset = calcLogmlChanges(inds, priorPar, data, cdat, logml, alpha, nCategories, nAttributes)

% Author(s): Paul Blomstedt, Pekka Marttinen

global COUNTS; 
global SUMCOUNTS;
global SAMPLEVAR; 
global SAMPLEMEAN; 
global PARTITION;   

npops = size(COUNTS,3);
muutokset = zeros(npops,1);
indsToBeMoved = inds;

if isempty(indsToBeMoved), return, end

i1 = PARTITION(indsToBeMoved(1));
diffInCounts = computeDiffInCounts(indsToBeMoved, size(COUNTS,2), data, nCategories);
diffInSumCounts = sum(diffInCounts);

COUNTS(:,:,i1) = COUNTS(:,:,i1) - diffInCounts;
SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:) - diffInSumCounts;
inds_i1 = find(PARTITION==i1); 
inds_i1_new = setdiff(inds_i1,indsToBeMoved); 
[SAMPLEMEAN(i1,:), SAMPLEVAR(i1,:)] = calcStats(inds_i1_new,cdat); 
updateLogmlTable(i1, priorPar, nCategories, nAttributes);

for i2 = 1:npops
    if i2 ~= i1
        COUNTS(:,:,i2) = COUNTS(:,:,i2) + diffInCounts;
        SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:) + diffInSumCounts;
        inds_i2 = find(PARTITION==i2);
        inds_i2_new = union(inds_i2,indsToBeMoved);
        [SAMPLEMEAN(i2,:), SAMPLEVAR(i2,:)] = calcStats(inds_i2_new,cdat); 
        updateLogmlTable(i2, priorPar, nCategories, nAttributes);
        logPosteriorNew = computeTotalLogPosterior(data,alpha);
        muutokset(i2) = logPosteriorNew - logml;
        COUNTS(:,:,i2) = COUNTS(:,:,i2) - diffInCounts;
        SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:) - diffInSumCounts;
        [SAMPLEMEAN(i2,:), SAMPLEVAR(i2,:)] = calcStats(inds_i2,cdat); 
        updateLogmlTable(i2, priorPar, nCategories, nAttributes);
    end
end
COUNTS(:,:,i1) = COUNTS(:,:,i1) + diffInCounts;
SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:) + diffInSumCounts;
[SAMPLEMEAN(i1,:), SAMPLEVAR(i1,:)] = calcStats(inds_i1,cdat); 
updateLogmlTable(i1, priorPar, nCategories, nAttributes);
