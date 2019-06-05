function [samplevar, samplemean] = ...
    initialStats(partition, cdat, npops, nAttributes)

% Author(s): Paul Blomstedt

samplevar = zeros(npops,nAttributes);
samplemean = zeros(npops,nAttributes);

for i=1:npops
    inds = find(partition==i);
    [samplemean(i,:), samplevar(i,:)] = calcStats(inds,cdat);
end

