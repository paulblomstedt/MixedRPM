function [partition, logPosterior]=...
    modelSearchGreedy(data,dist,npops, partition, roundTypes,discr,alpha)

% Author(s): Pekka Marttinen, Paul Blomstedt, Eero Linna

global PARTITION; global COUNTS;
global SUMCOUNTS; global LOGML_TABLE;
global SAMPLEVAR;
global SAMPLEMEAN;
global ADDITION_DIFFERENCE;
global REMOVAL_DIFFERENCE;
global JOIN_DIFFERENCE;
clearGlobalVars;
nObservations = size(data, 1);
nAttributes = size(data,2);
nCategories = length(discr)+1; % no. of discrete categories + one continuous category

% PRIOR VALUES:
priorPar = repmat(ones(nCategories, 1)./(nCategories), [1, nAttributes, npops]);

% Initialize PARTITION, COUNTS, SUMCOUNTS:
% Process data to continuous and discrete datasets so they can be dealt with separately.  
cdat = data;

% Calculate properties for discrete part
data = discretize(data, discr);
[sumcounts, counts] = initialCounts2(partition, data, npops, nCategories, nAttributes);                                                       
COUNTS = counts; SUMCOUNTS = sumcounts;

% Calculate properties for continuous part
cdat = discreteValuesTo0(cdat, discr);
cdat = normalizeData(cdat, nAttributes); % normalize data for each feature 
[samplevar, samplemean] = initialStats(partition, cdat, npops, nAttributes); 
SAMPLEVAR = samplevar; SAMPLEMEAN = samplemean;  
PARTITION = shiftdim(partition);  
clear partition; clear counts; clear sumcounts;
clear samplevar; clear samplemean; 


% Initialize LOGML_TABLE:
LOGML_TABLE = zeros(npops, 1);
updateLogmlTable(1:npops, priorPar, nCategories, nAttributes);

logPosterior = computeTotalLogPosterior(data, alpha);

REMOVAL_DIFFERENCE = zeros(nObservations,1);
REMOVAL_DIFFERENCE(:,:) = nan;
ADDITION_DIFFERENCE = zeros(nObservations,npops);
ADDITION_DIFFERENCE(:,:) = nan;
JOIN_DIFFERENCE = zeros(npops, npops);
JOIN_DIFFERENCE(:,:) = nan;

% ***********Doc:********************
% REMOVAL_DIFFERENCE(ind) tells the change in logml if ind is removed from
% its cluster. nan, if the cluster has changed, since the value was last
% calculated.
%
% ADDITION_DIFFERENCE(ind, pop) tells the change in logml if ind is added
% to cluster pop. nan, if the cluster has changed since the value was last
% calculated. Always nan, if pop is ind's own cluster.
%
% JOIN_DIFFERENCE(pop1,pop2) = tells the change in logml if pop1 and pop2
% are combined. nan, if either cluster has changed since the value was last
% calculated.
% ***********Doc end*****************

disp('Initialization:');
disp(['Partition: ' num2str(PARTITION')]);
disp(['Nclusters: ' num2str(length(unique(PARTITION)))]);
disp(['Log(ml*prior): ' num2str(logPosterior)]);
disp(' ');

% START SEARCH OF THE BEST PARTITION:

lever = zeros(1,14);
ready = 0;
while ready ~= 1
    for n = 1:length(roundTypes)
        round = roundTypes(n);
        magicnumber = 0;
        if  round==1 && lever(1)==0  % Moving an item into a different cluster
            inds = 1:nObservations;
            tmpTable = [inds' rand(nObservations,1)];
            tmpTable = sortrows(tmpTable,2);
            inds = tmpTable(:,1)';
            
            for ind = inds
                updateDifferenceTables(ind, priorPar, data, cdat, logPosterior, alpha, nCategories, nAttributes); %
                changes2 = REMOVAL_DIFFERENCE(ind) + ADDITION_DIFFERENCE(ind,:);
                changes2(PARTITION(ind)) = 0;
                [maxChange, index] = max(changes2,[],2);
                
                if maxChange>1e-5
                    updateGlobalVariables(ind, index, priorPar, data, cdat, nCategories, nAttributes); %
                    logPosterior = computeTotalLogPosterior(data, alpha);
                    
                    magicnumber = magicnumber+1;
                    lever = zeros(1,14);
                else
                    %disp(['Item ' num2str(ind) ' not moved.']);
                end
            end
            if magicnumber==0, lever(1)=1; end
            disp(['Step 1: ' num2str(magicnumber) ' individuals moved.']);
            
        elseif round==2 && lever(2)==0 % Merging two populations
            
            updateJoinDifference(priorPar, data, cdat, logPosterior, alpha, nCategories, nAttributes); %
            [maxChange, aux] = max(JOIN_DIFFERENCE(:));
            [i1, i2] = ind2sub([npops,npops],aux);
            
            if maxChange>1e-5
                changing = find(PARTITION==i1);
                if isempty(changing)
                    error('empty')
                end
                updateGlobalVariables(changing, i2, priorPar, data, cdat, nCategories, nAttributes);
                logPosterior = computeTotalLogPosterior(data, alpha);
                
                disp(['Step 2: Clusters ' num2str(i1) ' and ' num2str(i2) ' combined.']);
                lever = zeros(1,14);
            else
                disp('Step 2: no changes.');
                lever(2)=1;
            end
        elseif ismember(round, 3:4) && lever(round)==0  % Splitting a cluster
            nObservations = size(data,1);
            
            pops = 1:npops;
            tmpTable = [pops' rand(npops,1)];
            tmpTable = sortrows(tmpTable,2);
            pops = tmpTable(:,1)';
            
            splitpops = zeros(npops,1);
            for pop = pops
                
                maxChange = 0;
                inds2 = find(PARTITION==pop);
                nObservations2 = length(inds2);
                if nObservations2>4
                    
                    if round==3
                        dist3 = calcPartDist(inds2, dist, nObservations);
                        npops2 = min(20, floor(nObservations2 / 5));  % Split into npops2 parts
                    elseif round==4
                        dist3 = calcPartDist(inds2, dist, nObservations);
                        npops2 = 2;
                    end
                    
                    Z3 = linkage(dist3');
                    T3 = clusterInit(Z3, npops2);
                    
                    for i = 1:npops2
                        indsX = inds2(T3==i); indsX = indsX';
                        changes = calcLogmlChanges(indsX, priorPar, data, cdat, logPosterior, alpha, nCategories, nAttributes); %
                        [largest, index] = max(changes);
                        if largest>maxChange
                            maxChange = largest;
                            i2 = index;
                            changing = indsX;
                        end
                    end
                    if maxChange>1e-5
                        updateGlobalVariables(changing, i2, priorPar, data, cdat, nCategories, nAttributes); %
                        logPosterior = computeTotalLogPosterior(data, alpha);
                        splitpops(pop)=1;
                    end
                end
            end
            if ~isempty(find(splitpops, 1))
                disp(['Step ' num2str(round) ': ' num2str(length(find(splitpops))) ' populations were split.']);
                lever = zeros(1,14);
            else
                disp(['Step ' num2str(round) ': no changes.']);
                lever(round)=1;
            end
        end
    end
    
    roundTypes = [];
    
    if isempty(roundTypes)
        ready = 1;
    end
end

disp(' ');
disp('BEST PARTITION: ');
disp(['Nclusters: ' num2str(length(unique(PARTITION)))]);
disp(['Log(ml*prior): ' num2str(logPosterior)]);
disp(' ');

npops = removeEmptyClusters;  % Remove empty clusters

LOGML_TABLE = zeros(npops,1);
updateLogmlTable(1:npops, priorPar, nCategories, nAttributes);
logPosterior = computeTotalLogPosterior(data, alpha);

partition = PARTITION;

%----------------------------------------------------------------------------


function dist2 = calcPartDist(inds2, dist_orig, nObservations)
% Forms a subvector of the vector dist_orig, containing the distances
% between the items inds2
nObservations2 = length(inds2);
tmp = zeros(nchoosek(nObservations2,2),2);
rivi = 1;
for i=1:nObservations2-1
    for j=i+1:nObservations2
        tmp(rivi, 1) = inds2(i);
        tmp(rivi, 2) = inds2(j);
        rivi = rivi+1;
    end
end
tmp = (tmp(:,1)-1).*nObservations - tmp(:,1) ./ 2 .* (tmp(:,1)-1) + (tmp(:,2)-tmp(:,1));
dist2 = dist_orig(tmp);
