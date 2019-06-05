function [partition, logPosterior] = clusterMixed(data, nClustMax, discrVal, alpha)

% Performs Bayesian random partition clustering for a given data matrix of mixed discrete and continuous type.
%
% Input:
%
% data: data matrix, in which the rows correspond to the data items.
%
% nClustMax: estimated upper bound for the number of clusters. If the final 
% solution is close to the specified upper bound, the search should be repeated 
% with a higher upper bound.
%
% discrVal: vector of discrete values
%
% alpha (optional): concentration parameter for CRP prior, the default value is 1 
%
% Output:
%
% partition = estimated partition for the data items.
% logPosterior = log(marginal likelihood x prior)
%
% Author(s): Paul Blomstedt, Pekka Marttinen
%
% Reference:
% Paul Blomstedt, Jing Tang, Jie Xiong, Christian Granlund, and Jukka Corander. 
% A Bayesian predictive model for clustering data of mixed discrete and continuous type. 
% IEEE Transactions on Pattern Analysis and Machine Intelligence,
% 37(3):489-498, 2015.

if nargin < 4
    alpha = 1;
end

[data,discr] = mapDiscreteVals(data,discrVal);

[dist,Z] = calcEuclidDist(data);   % calculate distances
initPart = clusterInit(Z,nClustMax); % initial partition

% Use fixed sequence of search operators.
roundTypes = [2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 ...
            2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 ...
            2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 ...
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ...
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ...
            3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 ...
            3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 1 1 1 1 ...
            1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 ...
            1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 ...
            1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4];

% Start search from initial partition.
nClustMax = length(unique(initPart))+15;  % To be on the safe side, increase nClustMax.
partition = modelSearchGreedy(data,dist,nClustMax,initPart,roundTypes,discr,alpha);

% Repeat the search once more.
nClustMax = length(unique(partition))+5;
[partition, logPosterior] = modelSearchGreedy(data,dist,nClustMax,partition,roundTypes,discr,alpha);
