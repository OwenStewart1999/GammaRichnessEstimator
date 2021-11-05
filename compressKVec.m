function [kVec2, medK] = compressKVec(abundMat, kVec, logInd, nRuns, metric)
% compressKVec will take a kVec and abundMat and compress the kVec to a
% more reasonable set of values, thus reducing the computational load when
% determineKKMeans is called
% it will do so by applying  a given metric nRuns times, taking the median
% k value and creating a new kVec based on this value

% inputs:
% abundMat - a matrix of abundance vectors stored in column format - each
% column is an abundance vector for a different partition (i.e. quadrat)
% kVec - a vector contianing all k values to check through
% kVec - a vector contianing all k values to check through
% logInd - if equal to "log", abundMat will be log-transformed before
% clustering occurs using log(abundMat + 1)
% nRuns - number of times to run kMeans before taking the median
% metric - the specific metric used to determine each optimal k value
% collected - currently supports "Dunn's average" and "gap"

% output:
% kVec2 - the compressed kVec

% create a vector to store the different optimal k values in
kOptVec = zeros(1, nRuns);

% log transform abundMat if necessary
if logInd == "log"
    abundMat = log(abundMat + 1);
end

% determine the number of k values in the initial kVec
numK = length(kVec);

if metric == "Dunn's average"
    
    % now loop for nRuns - use Dunn's index using average distances, and do 3
    % replicates each for accuracy
    for r = 1:nRuns

        % calculate the pairwise distances between all points to make Dunn's
        % index more efficient
        D = pdist(abundMat');
        D = squareform(D);

        % create a vector to hold the Dunn's indices
        DunnIndVec = zeros(1, numK);

        % loop over each k value inside kVec, cluster, and calculate Dunn's
        % index
        for i = 1:numK

            % apply kmeans to the current number of clusters, with 3
            % replicates each, then calculate and store Dunn's index
            idx = kmeans(abundMat', kVec(i), 'Replicates', 3);
            DunnIndVec(i) = DunnsIndex(abundMat', idx, "average", [], D);

        end

        % choose the value of k which maximises Dunn's index - if kVec starts
        % at 1, remove this value as Dunn's index is not defined for one
        % cluster and gives an infinite value
        if kVec(1) == 1
            DunnIndVec(1) = 0;
        end
        [~, index] = max(DunnIndVec);
        kOptVec(r) = kVec(index);

    end

elseif metric == "gap"
    
    % remove k = 1 if present as this value is undefined for the gap
    % statistic
    if kVec(1) == 1
        kVec = kVec(2:end);
    end
    
    % if nRuns > 10 just set to 10 - 10 runs should be enough to get a
    % decent idea of the spread of values
    if nRuns > 10
        nRuns = 10;
    end
    
    % now loop for nRuns and apply the gap statistic at each run
    for r = 1:nRuns
        eval = evalclusters(abundMat', "kmeans", "gap", "KList", kVec, 'B', 50);
        kOptVec(r) = eval.OptimalK;
    end
    
end

% now find the median k value and ceil it in case it lies between 2 scores
medK = median(kOptVec);
medK = ceil(medK);

% now, create kVec2 based on medK
if medK >= 5
    
    % if the median k is >= 5, then go 1/2 k ciel'd above and below
    kVec2 = (medK - ceil(medK/2)):(medK + ceil(medK/2));
    
else
    
    % if the median k is < 5, then start at 2 and go to 5 + ceil(kMed/2)
    kVec2 = 2:(5 + ceil(medK/2));
    
end
end