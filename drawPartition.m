function drawPartition(part, data)

% Author(s): Pekka Marttinen, Paul Blomstedt

subplot(1,2,1)
image(data*95+5)
title('Original')

subplot(1,2,2)
nloci = size(data,2);

[X,I] = sort(part);
image(data(I,:)*95+5);

npops = length(unique(part));
last_inds = zeros(npops,1);

for i=1:npops
    last_inds(i) = find(X==i, 1, 'last' );
end

last_inds = last_inds+0.5;
line(repmat([0.5; nloci+0.5],[1 npops]),[last_inds';last_inds'],'Color','k'); 
title('Clustered')