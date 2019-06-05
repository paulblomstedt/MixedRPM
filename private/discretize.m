function data = discretize(data,discr)

% Author(s): Paul Blomstedt

ia =  ~ismember(data,discr);
data(ia) = max(discr)+1;