%% Demo of Bayesian random partition clustering for data of mixed discrete and continuous type.

% Reference:
% Paul Blomstedt, Jing Tang, Jie Xiong, Christian Granlund, and Jukka Corander. 
% A Bayesian predictive model for clustering data of mixed discrete and continuous type. 
% IEEE Transactions on Pattern Analysis and Machine Intelligence,
% 37(3):489-498, 2015.


% Load data
load examples

% Remark 1 
% The example data sets are for purposes of demonstration only, and are
% deliberately made easy for the clustering model to handle. They do not
% reflect the performance of the model for arbitrary real data.

% Remark 2
% Since the search algorithm is greedy, there is no guarantee of convergence 
% at the global maximum. It is therefore a good idea to run the algorithm multiple 
% times on a data set and monitor the log posterior value to see which solution 
% is the best one (i.e. has the highest value).

%% Example 1: 
% Very simple simulated zero-inflated dataset, which may still be difficult to 
% cluster correctly using a purely Gaussian or a purely binary clustering model.

% Cluster
partition = clusterMixed(dataSimple,10,0);

% Visualize
figdata = dataSimple./max(max(dataSimple)); % normalize for visual purposes
figure
drawPartition(partition,figdata)

% Evaluate similarity with known true partition
ARI = RandIndex(partSimple,partition);


%% Example 2
% Data with similar characteristics as real preprocessed amphetamine data.

% Cluster
partition = clusterMixed(dataReal,20,0);

% Visualize
figdata = log(dataReal+1)./max(max(log(dataReal+1))); % logarithmize and normalize for visual purposes
figure
drawPartition(partition,figdata)

%% Example 3
% Simulated data having two discrete categories.

% Cluster
partition = clusterMixed(data2, 20, [0 1]);
% Evaluate similarity with known true partition
ARI = RandIndex(part2,partition);

%% Example 4
% Simulated data having three discrete categories with arbitrary values.

% Cluster and evaluate (unnormalized) log posterior of a the obtained solution
[partition, logPosterior1] = clusterMixed(data3, 20, [5 23 42]);

% Evaluate similarity with known true partition
ARI1 = RandIndex(part3,partition);

% Evaluate adjusted Rand index and log posterior for an arbitrary partition
badPartition = datasample(1:10,size(data3,1)); 
ARI2 = RandIndex(part3,badPartition);
logPosterior2 = evaluate(data3,badPartition,[5 23 42]);


