function logPosterior = computeTotalLogPosterior(data, alpha)

% logPosterior = log(marginal likelihood) + log(prior)
% Author(s): Paul Blomstedt

global LOGML_TABLE;
global PARTITION;


ninds = size(data,1);
uniq = unique(PARTITION); 
n_c =histc(PARTITION,uniq);
k = numel(n_c); % Number of nonempty clusters
term = 0;
for c = 1:k;
    term = term + sum(log(1:(n_c(c)-1)));
end
logPosterior = sum(LOGML_TABLE,1)...
    + k*log(alpha) + term - sum(log(alpha + (1:ninds)-1));