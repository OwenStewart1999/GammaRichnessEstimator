function k = determineKKMeans(abundMat, kVec, logInd, metric, D)
% determineKKMeans will determine an optimal number of clusters (k) to use
% for k means clustering given the data in abundMat, the range of k values
% to trial in kVec and the clustering metric used in metric

% inputs:
% abundMat - a matrix of abundance vectors stored in column format - each
% column is an abundance vector for a different partition (i.e. quadrat)
% kVec - a vector contianing all k values to check through
% logInd - if equal to "log", abundMat will be log-transformed before
% clustering occurs using log(abundMat + 1)
% metric - specifies the metric used to determine the optimal k -
% "Silhouette" will use the maximum average silhouette scores between runs,
% "Dunn's average" will use Dunn's index using average distances and
% "Dunn's minmax" will use Dunn's index using maximum and minimum values
% "Dunn's average centroid" will use Dunn's index using average distances
% to the centroid - see the DunnsIndex function for more details
% "Dunn's average rep" will use the Dunn's average method, but repeat it 13 
% times and choose the highest estimate
% "AIC" or "BIC" will use the AIC or BIC assuming each quadrat's raw
% abundance data (without being log transformed) are the results given by a
% draw from one of k multinomial distributions, one for each cluster,
% where the mutlinomial probabilities are the average proportions of each
% species over the quadrats assigned to that cluster
% "gap", "Davies-Bouldin" and "Calinski-Harabasz" will use the gap
% statistic, Davies-Bouldin index and Calinski-Harabsz index respectively
% D - optional input, but would be a pairwise distance matrix in the same
% order as abundMat useful for increasing efficiency for Dunn's index
% methods

% output:
% k - the optimal number of clusters determined

% NOTE: this method assumes a correct input is passed into metric, and
% does not correct for incorrect metrics, as this is done inside the parent
% method gammaRichnessEstimator - if an unsupported metric is supplied,
% this method will simply return an error

% determine the number of k values to trial
numk = length(kVec);

% create a dummy empty D variable if not passed in
if nargin < 5 || isempty(D)
    D = [];
end

% determine which clustering metric to use
if strcmp(metric, "Silhouette")
    
    % if the silhouette metric is chosen, create a vector to hold the
    % silhouette scores in
    silhouetteVec = zeros(1, numk);
    
    % log transform the abundance data if necessary
    if strcmp(logInd, "log")
        abundMat = log(abundMat + 1);
    end
    
    % loop through the kVec and populate the silhouette score
    for i = 1:numk
        
        % apply kmeans to the current number of clusters, with 5
        % replicates each, then calculate and average the Silhouette score
        % and store in silhouetteVec
        idx = kmeans(abundMat', kVec(i), 'Replicates', 5);
        silhouetteVec(i) = mean(silhouette(abundMat', idx));
        
    end
    
    % choose the value of k which maximises the Silhouette score - if kVec
    % starts at 1, remove this value as the Silhouette Score is not defined
    % for one cluster and gives an infinite value
    if kVec(1) == 1
        silhouetteVec(1) = 0;
    end
    [~, index] = max(silhouetteVec);
    k = kVec(index);
    
elseif strcmp(metric, "Dunn's average")
    
    % if Dunn's index is chosen, create a vector to hold the Dunn's indices
    DunnIndVec = zeros(1, numk);
    
    % log transform the abundance data if necessary
    if strcmp(logInd, "log")
        abundMat = log(abundMat + 1);
    end
    
    % calculate the pairwise distances between each point once for
    % efficiency if D has not already been calculated
    if nargin < 5 || isempty(D)
        D = pdist(abundMat');
        D = squareform(D);
    end
    
    % loop through kVec and determine the Dunn's index for each point
    for i = 1:numk
        
        % apply kmeans to the current number of clusters, with 5
        % replicates each, then calculate and store Dunn's index
        idx = kmeans(abundMat', kVec(i), 'Replicates', 5);
        DunnIndVec(i) = DunnsIndex(abundMat', idx, "average", [], D);
        
    end
    
    % choose the value of k which maximises Dunn's index - if kVec starts
    % at 1, remove this value as Dunn's index is not defined for one
    % cluster and gives an infinite value
    if kVec(1) == 1
        DunnIndVec(1) = 0;
    end
    [~, index] = max(DunnIndVec);
    k = kVec(index);
    
elseif strcmp(metric, "Dunn's minmax")
    
    % if Dunn's index is chosen, create a vector to hold the Dunn's indices
    DunnIndVec = zeros(1, numk);
    
    % log transform the abundance data if necessary
    if strcmp(logInd, "log")
        abundMat = log(abundMat + 1);
    end
    
    % calculate the pairwise distances between each point once for
    % efficiency if D has not already been calculated
    if nargin < 5 || isempty(D)
        D = pdist(abundMat');
        D = squareform(D);
    end
    
    % loop through kVec and determine the Dunn's index for each point
    for i = 1:numk
        
        % apply kmeans to the current number of clusters, with 5
        % replicates each, then calculate and store Dunn's index
        idx = kmeans(abundMat', kVec(i), 'Replicates', 5);
        DunnIndVec(i) = DunnsIndex(abundMat', idx, "minmax", [], D);
        
    end
    
    % choose the value of k which maximises Dunn's index - if kVec starts
    % at 1, remove this value as Dunn's index is not defined for one
    % cluster and gives an infinite value
    if kVec(1) == 1
        DunnIndVec(1) = 0;
    end
    [~, index] = max(DunnIndVec);
    k = kVec(index);
    
elseif strcmp(metric, "Dunn's average centroid")
    
    % if Dunn's index is chosen, create a vector to hold the Dunn's indices
    DunnIndVec = zeros(1, numk);
    
    % log transform the abundance data if necessary
    if strcmp(logInd, "log")
        abundMat = log(abundMat + 1);
    end
    
    % calculate the pairwise distances between each point once for
    % efficiency if D has not already been calculated
    if nargin < 5 || isempty(D)
        D = pdist(abundMat');
        D = squareform(D);
    end
    
    % loop through kVec and determine the Dunn's index for each point
    for i = 1:numk
        
        % apply kmeans to the current number of clusters, with 5
        % replicates each, then calculate and store Dunn's index
        [idx, cent] = kmeans(abundMat', kVec(i), 'Replicates', 5);
        DunnIndVec(i) = DunnsIndex(abundMat', idx, "average centroid", cent, D);
        
    end
    
    % choose the value of k which maximises Dunn's index - if kVec starts
    % at 1, remove this value as Dunn's index is not defined for one
    % cluster and gives an infinite value
    if kVec(1) == 1
        DunnIndVec(1) = 0;
    end
    [~, index] = max(DunnIndVec);
    k = kVec(index);
    
elseif strcmp(metric, "AIC")
    
    % create a vector to hold AIC values
    AICVec = zeros(1, numk);
    
    % log transform abundance matrix if necessary
    if strcmp(logInd, "log")
        abundMatLog = log(abundMat + 1);
    end
    
    % loop through kVec and determine the AIC for each clustering scenario
    for i = 1:numk
        
        % apply kmeans to the current number of clusters, with 5
        % replicates each, then calculate and store the AIC
        if strcmp(logInd, "log")
            idx = kmeans(abundMatLog', kVec(i), 'Replicates', 5);
        else
            idx = kmeans(abundMat', kVec(i), 'Replicates', 5);
        end
        [AICVec(i), ~] = AICBICClust(abundMat, idx);
        
    end
    
    % choose the value of k which minimizes the AIC
    [~, index] = min(AICVec);
    k = kVec(index);
    
elseif strcmp(metric, "Dunn's average rep") 
    
    % if Dunn's index is chosen, create a vector to hold the Dunn's indices
    DunnIndVec = zeros(1, numk);
    
    % log transform the abundance data if necessary
    if strcmp(logInd, "log")
        abundMat = log(abundMat + 1);
    end
    
    % calculate the pairwise distances between each point once for
    % efficiency if D has not already been calculated
    if nargin < 5 || isempty(D)
        D = pdist(abundMat');
        D = squareform(D);
    end
    
    % this specific method involves repeating Dunn's Index using Average
    % distances and repeating it 13 times, so loop from 1:13
    maxK = 0;
    for rep = 1:13
    
        % loop through kVec and determine the Dunn's index for each point
        for i = 1:numk

            % apply kmeans to the current number of clusters, with only 2
            % replicates each, then calculate and store Dunn's index
            idx = kmeans(abundMat', kVec(i), 'Replicates', 2);
            DunnIndVec(i) = DunnsIndex(abundMat', idx, "average", [], D);

        end

        % choose the value of k which maximises Dunn's index - if kVec starts
        % at 1, remove this value as Dunn's index is not defined for one
        % cluster and gives an infinite value
        if kVec(1) == 1
            DunnIndVec(1) = 0;
        end
        [~, index] = max(DunnIndVec);
        kCurr = kVec(index);
        
        % check if the current k value found is greater than maxK, and if
        % so reassign maxK
        if kCurr > maxK
            maxK = kCurr;
        end
    
    end
    
    k = maxK;
    
elseif strcmp(metric, "BIC")
    
    % create a vector to hold BIC values
    BICVec = zeros(1, numk);
    
    % log transform abundance matrix if necessary
    if strcmp(logInd, "log")
        abundMatLog = log(abundMat + 1);
    end
    
    % loop through kVec and determine the BIC for each clustering scenario
    for i = 1:numk
        
        % apply kmeans to the current number of clusters, with 5
        % replicates each, then calculate and store the BIC
        if strcmp(logInd, "log")
            idx = kmeans(abundMatLog', kVec(i), 'Replicates', 5);
        else
            idx = kmeans(abundMat', kVec(i), 'Replicates', 5);
        end
        [~, BICVec(i)] = AICBICClust(abundMat, idx);
        
    end
    
    % choose the value of k which minimizes the BIC
    [~, index] = min(BICVec);
    k = kVec(index);
    
elseif strcmp(metric, "gap")
    
    % remove k = 1
    if kVec(1) == 1
        kVec = kVec(2:end);
    end
    
    % if "gap" is specified, apply the gap statistic using Matlab's inbuilt
    % implementation
    if logInd == "log"
        eval = evalclusters(log(abundMat' + 1), "kmeans", "gap", "KList", kVec, 'B', 50);
    else
        eval = evalclusters(abundMat', "kmeans", "gap", "KList", kVec);
    end
    k = eval.OptimalK;
    
elseif strcmp(metric, "Davies-Bouldin")
    
    % if "Davies-Bouldin" is specified, apply the Davies-Bouldin index
    % using Matlab's inbuilt implementation
    if logInd == "log"
        eval = evalclusters(log(abundMat' + 1), "kmeans", "DaviesBouldin", "KList", kVec);
    else
        eval = evalclusters(abundMat', "kmeans", "DaviesBouldin", "KList", kVec);
    end
    k = eval.OptimalK;
    
elseif strcmp(metric, "Calinski-Harabasz")
    
    % remove k = 1 from kVec if it appears
    if kVec(1) == 1
        kVec = kVec(2:end);
    end
    
    % if "Calinski-Harabasz" is specified, apply the Calinski-Harabasz
    % index using Matlab's inbuilt implementation
    if logInd == "log"
        eval = evalclusters(log(abundMat' + 1), "kmeans", "CalinskiHarabasz", "KList", kVec);
    else
        eval = evalclusters(abundMat', "kmeans", "CalinskiHarabasz", "KList", kVec);
    end
    k = eval.OptimalK;
    
else
    
    % return error if metric is not specified correctly
    error("Metric not correctly specified - please check inputs and read supported metrics")
    
end
end