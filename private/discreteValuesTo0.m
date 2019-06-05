function data = discreteValuesTo0(data,discr)

% Author(s): Paul Blomstedt

ia = ismember(data,discr);
data(ia) = 0;