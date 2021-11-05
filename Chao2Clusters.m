function alphaHatVec = Chao2Clusters(incMat, abundMat, idx)
% Chao2Cluster will calculate the estimated alpha richness of a set
% clusters using the Chao2 bias corrected method
% if the number of quadrats in a cluster is 2 or less, the Chao1 estimator
% will instead be used, as this is likely to be more accurate

% inputs:
% incMat - a matrix of incidence vectors stored in column format - each
% column is an incidence vector for a different partition (i.e. quadrat)
% abundMat - a matrix of abundance vectors stored in column format - each
% column is an abundance vector for a different partition (i.e. quadrat)
% idx - the index vector, which specifies the cluster each incidence vector
% belongs to

% output:
% alphaHatVec - a row vector containing the alpha richness estimators for
% each identified cluster/community - this will be stored in order of
% clusters identified in idx - i.e. the alpha hat estimator for community 5
% (as identified in idx) will be stored in alphaHatVec(5)

% first, determine the number of communities identified
k = max(idx);

% initialise the alphaHatVec
alphaHatVec = zeros(1, k);

% loop through each community
for i = 1:k
    
    % create a mask indicating which incidence vectors belong to community
    % i
    mask = (idx == i);
    
    % if the number of quadrats in the cluster is 2 or less, apply the
    % Chao1 method instead
    if sum(mask) <= 2
        
        % calculate the Chao1 estimated species richness
        alphaHatVec(i) = Chao1BiasCorrected(sum(abundMat(:, mask), 2));
        
    else
    
        % calculate the Chao2 estimated species richness using the
        % incidence vectors associated with the current cluster/community
        alphaHatVec(i) = Chao2BiasCorrected(incMat(:, mask));
    
    end
    
end

end

