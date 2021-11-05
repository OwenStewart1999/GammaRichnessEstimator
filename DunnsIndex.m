function DI = DunnsIndex(X, idx, metric, cent, D)
% DunnsIndex will calculate the Dunn's Index for cluster performance, given
% observations in X and a clustering specified by idx
% this method can use a number of different versions of Dunn's Index, which
% can be specified by the metric input

% inputs:
% X - the set of observations to be applied to, where each row is a
% different observation, and each column is a different variable
% idx - the index vector, which specifies the cluster each observation
% belongs to
% metric - the metric used to calculate Dunn's Index - currently three
% options are supported: "average", "minmax", "average centroid"
% "average" uses the average pairwise distances between points in either
% cluster for the inter - cluster distance, and the average pairwise
% distance between all points in a cluster for the intra - cluster distance
% "minmax" uses the minimum pairwise distance for the inter - cluster
% distance, and uses the maximum pairwise distance for the intra - cluster
% distance
% "average centroid" uses the average pairwise distances between points in
% either cluster for the inter - cluster distance, and the average distance
% to the centroid from each point for the intra - cluster distance
% "average efficient" will be almost identical to to "average", but will
% use some shortcuts to reduce the computational strain (doesn't actually
% help much though, so don't bother using this version)
% cent - a matrix containing the cluster centroids for each cluster, which
% can be outputted from the kmeans method itself, each row corresponds to a
% different centroid however is an optional input as it can be calculated
% later instead
% D - optional input, but can be used to save efficiency - represents a
% pairwise distance matrix produced by pdist() where rows and columns are
% in the same order as X

% output:
% DI - the Dunn's Index for this particular clustering arrangement

% determine the number of clusters
k = length(unique(idx));

% determine the number of observations in each cluster
nObsVec = zeros(1, k);
for i = 1:k
    nObsVec(i) = sum(idx == i);
end

% re - order the observations, so all observations from cluster 1 are
% first, all observations from cluster 2 are next etc
[~, I] = sort(idx);
X = X(I, :);

% if D has not been passed in as an argument, calculate D
if nargin < 5 || isempty(D)
    % determine the pairwise distances between each point, then convert to
    % a square matrix
    D = pdist(X);
    D = squareform(D);
else
    % otherwise if D has already been passed in, need to rearrange so that
    % it is in the correct order based on the new indexing
    D = D(I, I);
end

% create a vector detailing where each cluster starts, and another where
% each cluster ends
startInd = zeros(1, k);
startInd(1) = 1;
for i = 2:k
    startInd(i) = startInd(i-1) + nObsVec(i-1);
end
endInd = startInd + nObsVec - 1;

% use if statements to alter the metric used to determine distances
if metric == "average"
    
    % determine the maximum intra-cluster distance using the mean of all
    % pairwise distances in each cluster
    maxIntraCDist = 0;
    for i = 1:k
        if nObsVec(i) > 1
            % if the number of observations is greater than one, determine
            % the average pairwise distance between all points in the
            % current cluster
            intraCDistCurr = sum(sum(D(startInd(i):endInd(i), startInd(i):endInd(i))));
            intraCDistCurr = intraCDistCurr / (nObsVec(i)^2 - nObsVec(i));
            if intraCDistCurr > maxIntraCDist
                maxIntraCDist = intraCDistCurr;
            end
        end
    end

    % determine the minimum inter-cluster distance, using the average
    % pairwise distance between all points in each cluster - this value
    % will be determined for all possible pairs of clusters
    minInterCDist = inf;
    for i = 1:(k - 1)
        for j = (i + 1):k

            % determine the average pairwise distance between all points in
            % cluster i and cluster j
            interCDistCurr = sum(sum(D(startInd(i):endInd(i), startInd(j):endInd(j))));
            interCDistCurr = interCDistCurr / (nObsVec(i) * nObsVec(j));

            % check if the current inter-cluster distance is the minimum
            if interCDistCurr < minInterCDist
                minInterCDist = interCDistCurr;
            end
        end
    end

elseif metric == "minmax"
    
    % determine the maximum intra-cluster distance minimum of each
    % cluster's pairwise distances
    maxIntraCDist = 0;
    for i = 1:k
        if nObsVec(i) > 1
            % if the number of observations is greater than one, determine
            % the max pairwise distance between all points in the
            % current cluster
            intraCDistCurr = max(max(D(startInd(i):endInd(i), startInd(i):endInd(i))));
            if intraCDistCurr > maxIntraCDist
                maxIntraCDist = intraCDistCurr;
            end
        end
    end

    % determine the minimum inter-cluster distance, using the minimum
    % pairwise distance between all points in each cluster - this value
    % will be determined for all possible pairs of clusters
    minInterCDist = inf;
    for i = 1:(k - 1)
        for j = (i + 1):k

            % determine the minimum pairwise distance between all points in
            % cluster i and cluster j
            interCDistCurr = min(min(D(startInd(i):endInd(i), startInd(j):endInd(j))));

            % check if the current inter-cluster distance is the minimum
            if interCDistCurr < minInterCDist
                minInterCDist = interCDistCurr;
            end
        end
    end
    
elseif metric == "average centroid"
    
    % check if cent has been supplied by the user, and if not, calculate
    % the centroids of each cluster
    if nargin < 4 || isempty(cent)
        
        % initialise cent, where each row is a different cluster's centroid
        dim = size(X, 2);
        cent = zeros(k, dim);
        for i = 1:k
            cent(i, :) = sum(X(startInd(i):endInd(i), :), 1) / nObsVec(i);
        end
        
    end
    
    % determine the maximum intra-cluster distance using the distance to
    % each cluster's centroid
    maxIntraCDist = 0;
    for i = 1:k
        
        % calculate the average distance from the cluster centroid to each
        % point in the cluster
        intraCDistCurr = sum(pdist2(cent(i, :), X(startInd(i):endInd(i), :))) / nObsVec(i);
        if intraCDistCurr > maxIntraCDist
            maxIntraCDist = intraCDistCurr;
        end
    end
    
    % determine the minimum inter-cluster distance, using the average
    % pairwise distance between all points in each cluster - this value
    % will be determined for all possible pairs of clusters
    minInterCDist = inf;
    for i = 1:(k - 1)
        for j = (i + 1):k

            % determine the average pairwise distance between all points in
            % cluster i and cluster j
            interCDistCurr = sum(sum(D(startInd(i):endInd(i), startInd(j):endInd(j))));
            interCDistCurr = interCDistCurr / (nObsVec(i) * nObsVec(j));

            % check if the current inter-cluster distance is the minimum
            if interCDistCurr < minInterCDist
                minInterCDist = interCDistCurr;
            end
        end
    end
    
elseif metric == "average efficient"
    
    % for some reason this isn't even more efficient so honestly don't
    % bother using it
    
    % check if cent has been supplied by the user, and if not, calculate
    % the centroids of each cluster
    if nargin < 4 || isempty(cent)
        
        % initialise cent, where each row is a different cluster's centroid
        dim = size(X, 2);
        cent = zeros(k, dim);
        for i = 1:k
            cent(i, :) = sum(X(startInd(i):endInd(i), :), 1) / nObsVec(i);
        end
        
    end

    % determine the maximum intra-cluster distance using the mean of all
    % pairwise distances in each cluster
    maxIntraCDist = 0;
    for i = 1:k
        if nObsVec(i) > 1
            % if the number of observations is greater than one, determine
            % the average pairwise distance between all points in the
            % current cluster
            intraCDistCurr = sum(pdist(X(startInd(i):endInd(i), :))) / (nchoosek(nObsVec(i), 2));
            if intraCDistCurr > maxIntraCDist
                maxIntraCDist = intraCDistCurr;
            end
        end
    end
    
    % determine the pairwise distance between all cluster centroids
    D2 = pdist(cent);
    
    % figure out how many community pairs to investigate - investigate the
    % bottom 20% if k >= 9, if k < 9 choose from specified values
    if k < 9
        nInvVec = [0 1 3 4 6 7 7 7];
        nInv = nInvVec(k);
    else
        nInv = ceil(0.2*nchoosek(k, 2));
    end
    
    % find the nInv closest pairs of clusters
    [~, I2] = mink(D2, nInv);
    
    % need to now convert I2 into the actual communities it represents - do
    % this by creating 2 arrays which show the first and second community
    % in the distance pairing
    C1Ind = [];
    for i = 1:(k-1)
        C1Ind = [C1Ind i*ones(1, k-i)];
    end
    C2Ind = [];
    for i = 1:(k-1)
        C2Ind = [C2Ind (1 + i):k];
    end
    
    % determine the minimum inter-cluster distance, using the average
    % pairwise distance between all points in each cluster, and only do so
    % for the closest nInv clusters
    minInterCDist = inf;
    for i = 1:nInv
        
        % determine the average pairwise distance between all points in
        % current 2 clusters being examined
        C1 = C1Ind(I2(i));
        C2 = C2Ind(I2(i));
        interCDistCurr = sum(sum(pdist2(X(startInd(C1):endInd(C1), :), X(startInd(C2):endInd(C2), :)))) / (nObsVec(C1) * nObsVec(C2));

        % check if the current inter-cluster distance is the minimum
        if interCDistCurr < minInterCDist
            minInterCDist = interCDistCurr;
        end
        
    end
    
else
    % if metric passed in does not match any pre - determined options, set
    % the DI to 0
    minInterCDist = 0;
    maxIntraCDist = 1;
end

% calculate Dunn's Index
DI = minInterCDist / maxIntraCDist;

% if for whatever reason, the DI comes as inf just set it to 0 instead, as
% the DI should be maximised
if DI == inf
    DI = 0;
end

end