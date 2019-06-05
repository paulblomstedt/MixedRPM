function [mappedData,newDiscr] = mapDiscreteVals(data,discr)

% Author(s): Paul Blomstedt

mappedData = data;
discr = unique(discr);
for i = 1:length(discr)
    mappedData(data==discr(i)) = i-1;
end
newDiscr = 0:length(discr)-1;