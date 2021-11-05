function propMat = calcProp(abundMat, idx)
% calcProp will calculate the proportion vectors for each of the
% communities, which will be used in the likelihood calculation for a given
% clustering
% proportion vectors are simply the proportions of each species seen in
% each community

% inputs:
% abundMat - a matrix of abundance vectors stored in column format - each
% column is an abundance vector for a different partition (i.e. quadrat)
% idx - the index vector, specifying which community each abundance vector
% (i.e. quadrat) belongs to
% this method assumes that idx begins at 1 and increases through every
% integer until k, the number of communities

% outputs:
% propMat - a matrix holding the underlying proportion vectors for
% each community, where each column refers to a different community (in the
% same order as idx)

% determine the number of species, communities and quadrats
nSpec = size(abundMat, 1);
k = length(unique(idx));

% initialise propMat
propMat = zeros(nSpec, k);

% loop through each community and create the proportion vectors for each
% community
for i = 1:k
    
    % determine which quadrats belong to the current community
    mask = (idx == i);
    
    % sum abundance vectors belonging to this community
    sumVec = sum(abundMat(:, mask), 2);
    
    % convert to proportions and store in propMat
    propMat(:, i) = sumVec / sum(sumVec);
    
end

end