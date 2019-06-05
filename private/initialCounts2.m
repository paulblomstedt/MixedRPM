function [sumcounts, counts] = ...
    initialCounts2(partition, mappedData, npops, nCategories, nAttributes)

% Initializes counts and sumcounts

% Author(s): Paul Blomstedt, Pekka Marttinen

counts = zeros(nCategories,nAttributes,npops);
sumcounts = zeros(npops,nAttributes);

for i=1:npops
    
    inds = find(partition==i);
    pop_size = length(inds); 
    sumcounts(i,:) = pop_size;
        
    for j=1:nAttributes
        for k=1:nCategories
            counts(k, j, i) = sum(mappedData(inds,j)==(k-1));
        end
    end
end


