function [AIC, BIC] = AICBICClust(abundMat, idx)
% AICBICClust will calculate the AIC and BIC for a given clustering
% assignment, assuming each cluster is a set of quadrats, each of which has
% abundances drawn from a latent multinomial distribution specific to that
% cluster
% this assumes that the probabilities for each multinomial event (i.e. a
% species appears) are the proportions of those species' abundances over
% the entire community

% inputs:
% abundMat - a matrix of abundance vectors stored in column format - each
% column is an abundance vector for a different partition (i.e. quadrat)
% idx - the index vector, which specifies the cluster each abundance vector
% belongs to

% outputs:
% AIC - the Akaike Information Criterion
% BIC - the Bayesian Information Criterion

% first, calculate the log-likelihood of the given cluster
LL = logLikelihoodCluster(abundMat, idx);

% determine the number of clusters, observations, and dimensions of the
% data
k = max(idx);
obs = size(abundMat, 2);
dim = size(abundMat, 1);

% determine the AIC - assume the number of estimated parameters is the
% number of clusters multiplied by the dimensions in the data (i.e.
% species)
AIC = 2*k*dim - 2*LL;

% determine the BIC - assume again that the number of estimated parameters
% is the number of clusters multiplied by the dimensions in the data
BIC = k*dim*log(obs) - 2*LL;

end