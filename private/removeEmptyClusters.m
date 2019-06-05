function npops = removeEmptyClusters
% Removes empty populations from COUNTS and SUMCOUNTS 
% Updates npops and PARTITION.

% Author(s): Paul Blomstedt, Pekka Marttinen

global COUNTS;
global SUMCOUNTS;
global PARTITION;

notEmpty = find(any(SUMCOUNTS,2));
COUNTS = COUNTS(:,:,notEmpty);
[x, y, z] = size(COUNTS);
COUNTS(:,:,z+1) = zeros(x,y); 
SUMCOUNTS = SUMCOUNTS(notEmpty,:);
SUMCOUNTS(z+1,:) = zeros(1,y); 

for n=1:length(notEmpty)
    tmp = PARTITION==notEmpty(n);
    PARTITION(tmp)=n;
end
npops = length(notEmpty)+1; 

