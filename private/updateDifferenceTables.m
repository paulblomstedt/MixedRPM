function updateDifferenceTables(ind, priorPar, data, cdat, logml, alpha, nCategories, nAtttributes)

% Author(s): Paul Blomstedt, Pekka Marttinen

global COUNTS;      
global SUMCOUNTS;
global PARTITION;   
global SAMPLEVAR; 
global SAMPLEMEAN; 
global ADDITION_DIFFERENCE;
global REMOVAL_DIFFERENCE;

rem_old = REMOVAL_DIFFERENCE;
add_old = ADDITION_DIFFERENCE;


diffInCounts = computeDiffInCounts(ind, size(COUNTS,2), data, nCategories);
diffInSumCounts = sum(diffInCounts);

if isnan(rem_old(ind))
    % Update removal difference for the individual:
    i1 = PARTITION(ind);
    
    COUNTS(:,:,i1) = COUNTS(:,:,i1) - diffInCounts;
    SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:) - diffInSumCounts;
    inds_i1 = find(PARTITION==i1); 
    inds_i1_new = setdiff(inds_i1,ind); 
    [SAMPLEMEAN(i1,:), SAMPLEVAR(i1,:)] = calcStats(inds_i1_new,cdat); 
    updateLogmlTable(i1, priorPar, nCategories, nAtttributes);
    logPosteriorNew = computeTotalLogPosterior(data,alpha);
    rem_old(ind) = logPosteriorNew-logml;
    COUNTS(:,:,i1) = COUNTS(:,:,i1) + diffInCounts;
    SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:) + diffInSumCounts;
    [SAMPLEMEAN(i1,:), SAMPLEVAR(i1,:)] = calcStats(inds_i1,cdat);
    updateLogmlTable(i1, priorPar, nCategories, nAtttributes);
end

new_pops = isnan(add_old(ind,:));
new_pops(PARTITION(ind)) = 0;   
new_pops = find(new_pops);

for i2 = new_pops
    % Update addition differences for the individual:
    COUNTS(:,:,i2) = COUNTS(:,:,i2) + diffInCounts;
    SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:) + diffInSumCounts;
    inds_i2 = find(PARTITION==i2); 
    inds_i2_new = union(inds_i2,ind);
    [SAMPLEMEAN(i2,:), SAMPLEVAR(i2,:)] = calcStats(inds_i2_new,cdat); 
    updateLogmlTable(i2, priorPar, nCategories, nAtttributes);
    logPosteriorNew = computeTotalLogPosterior(data,alpha);
    add_old(ind,i2) = logPosteriorNew - logml;
    COUNTS(:,:,i2) = COUNTS(:,:,i2) - diffInCounts;
    SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:) - diffInSumCounts;
    [SAMPLEMEAN(i2,:), SAMPLEVAR(i2,:)] = calcStats(inds_i2,cdat); 
    updateLogmlTable(i2, priorPar, nCategories, nAtttributes);
end

REMOVAL_DIFFERENCE = rem_old;
ADDITION_DIFFERENCE = add_old;
