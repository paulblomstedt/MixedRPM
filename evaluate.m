function logPosterior = evaluate(data,partition,discrVal)
% Calculates the marginal likelihood for the data using a given partition. 

% Author(s): Pekka Marttinen, Paul Blomstedt

global PARTITION; 
global COUNTS;
global SUMCOUNTS; 
global LOGML_TABLE;
global SAMPLEVAR;
global SAMPLEMEAN;
clearGlobalVars;

if nargin < 4
    alpha = 1;
end

[data,discr] = mapDiscreteVals(data,discrVal);

nAttributes = size(data,2);
nCategories = length(discr)+1; % no. of discrete categories + one continuous category
npops = length(unique(partition));

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
% LOGML_TABLE = zeros(npops, 1);
% updateLogmlTable(1:npops, priorPar, nCategories, nAttributes);
% 
% logPosterior = computeTotalLogPosterior(data, alpha);


LOGML_TABLE = zeros(npops,1);
updateLogmlTable(1:npops, priorPar, nCategories, nAttributes);
logPosterior = computeTotalLogPosterior(data, alpha);