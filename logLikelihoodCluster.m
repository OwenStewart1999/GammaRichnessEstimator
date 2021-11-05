function LL = logLikelihoodCluster(abundMat, idx)
% logLikelihoodCluster will calculate the log-likelihood of a clustering
% configuration given the set of abundance vectors and indices
% this assumes that each cluster is a set of observations taken from a
% different multinomial distribution

% inputs:
% abundMat - a matrix of abundance vectors stored in column format - each
% column is an abundance vector for a different partition (i.e. quadrat)
% idx - the index vector, which specifies the cluster each abundance vector
% belongs to

% outputs:
% LL - the log-likelihood

% determine the number of clusters
k = max(idx);

% calculate the proportion/probability vectors for each community
propMat = calcProp(abundMat, idx);

% add up the log-likelihoods from each of the communities
LL = 0;
for i = 1:k
    
    % determine the vectors belonging to this community, then calculate the
    % log-likelihood of that community
    mask = idx == i;    
    LL = LL + logLikelihoodSingleComm(abundMat(:, mask), propMat(:, i));
    
end
end