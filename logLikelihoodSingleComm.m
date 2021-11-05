function LL = logLikelihoodSingleComm(abundMat, probVec)
% logLikelihoodSingleComm will calculate the log likelihood of observing a
% set of abundance vectors stored in abundMat assuming each abundance
% vector is a draw from a multinomial distribution with probabilities
% stored in probVec 
% this is under the assumption probVec is the average proportions of each
% species in the community

% inputs:
% abundMat - a matrix of abundance vectors stored in column format - each
% column is an abundance vector for a different partition (i.e. quadrat)
% probVec - a vector describing the average proportions of each species in
% the current community, i.e. probVec(i) is the proportion of individuals
% in the community of species i (equivalent to the probability of
% choosing that species from in a random sample of the population)

% outputs:
% LL - the log-likelihood value calculated

% to avoid later issues with log(0), remove any species which have a 0
% probability of appearing in the community
mask = probVec == 0;
abundMat = abundMat(~mask, :);
probVec = probVec(~mask);

% determine the number of species with a nonzero probability, and the
% number of quadrats in this community
s = length(probVec);
N = size(abundMat, 2);

% create a vector of the n values (i.e. number of individuals in each
% quadrat)
nVec = sum(abundMat, 1);

% initialise log-likelihood
LL = 0;

% first, account for the outcomes multiplied by the log of their
% probabilities
for j = 1:s
    
    LL = LL + sum(abundMat(j, :)) * log(probVec(j));
        
end

% second, account for the last 2 summation terms - the log of n! and the
% sum of the logs of each observation
for i = 1:N
    
    % first determine if the n value for the current quadrat exceeds 170,
    % where matlab encounters overflow and treat cases differently
    if nVec(i) > 170
        % use log laws to convert the log of a factorial into a sum of logs
        for y = 2:nVec(i)
            LL = LL + log(y);
        end
    else
        % otherwise can just use factorial
        LL = LL + log(factorial(nVec(i)));
    end
    
    % now subtract away the log of the factorial of each outcome
    for j = 1:s
        % check again if numerical overflow will occur
        if abundMat(j, i) > 170
            for y = 2:abundMat(j, i)
                LL = LL - log(y);
            end
        else
            LL = LL - log(factorial(abundMat(j, i)));
        end
    end
end

end