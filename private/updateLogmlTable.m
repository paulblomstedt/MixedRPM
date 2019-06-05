function updateLogmlTable(pops, priorPar, nCategories, nAttributes)

% Updates global LOGML_TABLE variable. 
% LOGML_TABLE is an npops*1 array, which includes logml-values for different
% populations. This function updates the logml values for given populations
% "pops", based on the observed counts for the populations, given in global
% variables COUNTS, SUMCOUNTS, SAMPLEVAR and SAMPLEMEAN.

% Author(s): Paul Blomstedt, Pekka Marttinen

global COUNTS; 
global SUMCOUNTS;
global SAMPLEVAR; 
global SAMPLEMEAN; 
global LOGML_TABLE;
z = length(pops);

% discrete part
counts = COUNTS;
sumcounts = SUMCOUNTS;

sumPrior = sum(priorPar(:,:,pops),1); 
if z>1 
    sumPrior = squeeze(sumPrior); 
    sumPrior = sumPrior'; 
end 

term = gammaln(counts(:,:,pops)+priorPar(:,:,pops));  
term = reshape(term, [nCategories, nAttributes, z]); 
term = sum(term,1); 
if z>1 
    term = squeeze(term); 
    term = term'; 
end

sumGammaPrior = sum(gammaln(priorPar(:,:,pops)));
if z>1 
    sumGammaPrior = squeeze(sumGammaPrior); 
    sumGammaPrior = sumGammaPrior'; 
end

discr = gammaln(sumPrior) - gammaln(sumPrior+sumcounts(pops,:)) + term - sumGammaPrior;

% continuous part
samplevar = SAMPLEVAR;
samplemean = SAMPLEMEAN;

n_2cj = counts(nCategories,:,pops);
if z>1 
    n_2cj = squeeze(n_2cj); 
    n_2cj = n_2cj'; 
end
alpha_2cj = 1+n_2cj./2;
beta_2cj = 1+n_2cj./2.*samplevar(pops,:)+n_2cj.*samplemean(pops,:).^2./(2.*(1+n_2cj));
rho_2cj = 1+n_2cj;

cont = gammaln(alpha_2cj)-alpha_2cj.*log(beta_2cj)...
    -1/2.*log(rho_2cj)-n_2cj./2.*log(2*pi);
cont(n_2cj==0)=0; % zero data should not contribute to logml

% log(marginal likelihood)

table = discr+cont;
table = sum(table,2); % sum over dimensions
LOGML_TABLE(pops) = table;










