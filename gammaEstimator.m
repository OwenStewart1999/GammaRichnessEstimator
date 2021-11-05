function gamma = gammaEstimator(alphaHatVec, bCell)
% gammaEstimator will provide the species estimate for the entire region by
% exploiting the underlying community structure
% this method only requires the alpha richness for each community to be
% already computed, and does not require this to be done a certain way -
% hence alpha's can be estimated using any different method

% inputs:
% alphaHatVec - a row vector containing the alpha richness estimators for
% each identified cluster/community
% bCell - a cell array which stores incidence vectors for each community in
% an indexed format - i.e. if a community contains species 2, 3, and 5 out
% of 5 known species, bCell will hold the vector [2 3 5] for that
% community, rather than [0 1 1 0 1]
% bCell will be ordered community by community, i.e. bCell{3} will
% provide the indexed incidence vector for community 3
% the indexed incidence vectors will be stored as column vectors, the only
% time bCell is used is with the intersect command in bFunc, for which the
% orientation actually doesn't matter anyway
% bCell MUST have the communities represented in the same order as
% alphaHatVec, i.e. alphaHatVec(3) and bCell{3} should both refer to the
% respective alpha estimate and indexed incidence vector for community 3

% output:
% gamma - the estimated species richness for the entire region

% determine the number of communities
N = length(alphaHatVec);

% initialise gamma
gamma = 0;

% calculate gamma using the equation listed in the manuscript for gamma
% richness - note, that bFunc must be divided by alpha_i as bFunc only
% calculates the intersections and their lengths, but does not scale by
% alpha_i
for i = 1:(N - 1)
    
    % calculate the sum result of all beta's necessary for community i
    betaSum = bFunc(bCell{i}, i, 2, bCell, N);
    gamma = gamma + ...
        (1 +  betaSum/alphaHatVec(i)) * alphaHatVec(i);
end

% add the alpha estimate for the final community
gamma = gamma + alphaHatVec(N);

end