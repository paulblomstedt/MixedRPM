function updateJoinDifference(priorPar, data, cdat, logml,alpha, nCategories, nAttributes) 

% Author(s): Paul Blomstedt, Pekka Marttinen

global COUNTS;      global SUMCOUNTS;
global PARTITION;   global LOGML_TABLE;
global JOIN_DIFFERENCE;
global SAMPLEVAR; 
global SAMPLEMEAN; 
npops = size(COUNTS,3);

for i1 = 1:npops-1
    indsToBeMoved = find(PARTITION==i1);
    if isempty(indsToBeMoved)
        % Cluster i1 is empty
        JOIN_DIFFERENCE(i1,i1+1:npops) = 0;
        JOIN_DIFFERENCE(i1+1:npops,i1) = 0;
    else
        diffInCounts = computeDiffInCounts(indsToBeMoved, size(COUNTS,2), data, nCategories);
        diffInSumCounts = sum(diffInCounts);
        
        unknown_pops = find(isnan(JOIN_DIFFERENCE(i1,i1+1:end)));
        unknown_pops = unknown_pops+i1;
        
        LOGML_TABLE(i1) = 0; % the result of removing all individuals from cluster i1 
        
        for i2 = unknown_pops
            COUNTS(:,:,i2) = COUNTS(:,:,i2) + diffInCounts;
            SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:) + diffInSumCounts;
            inds_i2 = find(PARTITION==i2); 
            inds_i2_new = union(inds_i2,indsToBeMoved); 
            [SAMPLEMEAN(i2,:), SAMPLEVAR(i2,:)] = calcStats(inds_i2_new,cdat); 
            updateLogmlTable(i2, priorPar, nCategories, nAttributes);
            logPosteriorNew = computeTotalLogPosterior(data, alpha);
            JOIN_DIFFERENCE(i1,i2) = logPosteriorNew-logml;
            JOIN_DIFFERENCE(i2,i1) = logPosteriorNew-logml;
            COUNTS(:,:,i2) = COUNTS(:,:,i2) - diffInCounts;
            SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:) - diffInSumCounts;
            [SAMPLEMEAN(i2,:), SAMPLEVAR(i2,:)] = calcStats(inds_i2,cdat);
            updateLogmlTable(i2, priorPar, nCategories, nAttributes);
        end
        updateLogmlTable(i1, priorPar, nCategories, nAttributes);
    end
end