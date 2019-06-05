function [localMean, localVar] = calcStats(inds,cdat)
% CALCSTATS Calculates local mean value and local variance for the given
% indexes of cdat. Ignores zero values. All discrete values are
% mapped to zero at the beginning of model_search_parallel so
% CALCSTATS ignores all discrete values.
%
% Usage: [localMean, localVar] = CALCSTATS(inds, cdat)
%
% Input values: 
% inds is a row vector containing the indexes of cdat for which mean value and variance are
% to be calculated
% cdat is the data for which local mean and variance are calculated.
%
% Output values:
% localMean is a vector of mean values of each column

% Variables used in the code: 
%
% nContinuous is a vector that tells how many variables in each column of cdat in rows specified by inds are continuos.

% Author(s): Paul Blomstedt, Eero Linna

nContinuous = sum(cdat(inds,:)~=0);

% Avoid dividing by zero. 
zero = nContinuous==0;
nContinuous(zero) = 1;
% If the i:th value in nc is zero it means all the variables in column i in rows specified by inds were discrete

sumData = cdat(inds,:);
localMean = sum(sumData)./nContinuous;
localVar = sum(sumData.^2)./nContinuous-localMean.^2;




