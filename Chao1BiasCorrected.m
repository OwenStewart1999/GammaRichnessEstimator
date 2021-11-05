function Chao1 = Chao1BiasCorrected(abundVec)
% Chao1BiasCorrected will calculate the bias corrected version of the Chao1
% estimator given a vector of species abundances
% this method assumes that the vector may still have 0's in it, which could
% occur as the result of a vector being taken from a sample of a larger
% area, in which there are a known number of species

% input:
% abundMat - a matrix of abundance vectors stored in column format - each
% column is an abundance vector for a different partition (i.e. quadrat)

% output:
% Chao1 - the Chao1 estimated species richness

% first, determine the number of singletons and doubletons
f1 = sum(abundVec == 1);
f2 = sum(abundVec == 2);
Sobs = sum(abundVec > 0);

% calculate the bias corrected Chao1 estimator
Chao1 = Sobs + (f1 * (f1 - 1)) / (2 * (f2 + 1));

end

