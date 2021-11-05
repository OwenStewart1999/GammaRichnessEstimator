function Chao2 = Chao2BiasCorrected(incMat)
% Chao2BiasCorrected will calculate the bais corrected Chao2 estimate for
% alpha diversity, using incMat, which will store the incidence vectors
% required

% input:
% incMat - a matrix of incidence vectors stored in column format - each
% column is an incidence vector for a different partition (i.e. quadrat)

% output:
% Chao2 - the estimate for alpha richness using the Chao2 estimator in its
% bias corrected form

% determine then number of samples/quadrats
m = size(incMat, 2);

% sum along rows to determine the number of incidences of each species
repIncVec = sum(incMat, 2);

% determine the total number of species observed in the collection of
% quadrats, the total number of unique species and the total number of
% duplicate species
Sobs = sum(repIncVec > 0);
q1 = sum(repIncVec == 1);
q2 = sum(repIncVec == 2);

% calculate the Chao2 estimate for estimated species richness
Chao2 = Sobs + ((m - 1) * q1 * (q1 - 1)) / (m * 2 * (q2 + 1));

end