function [estimateVec, optKVec] = gammaRichnessEstimator(abundMat, nRuns, logInd, kVec, replicates, metric)
% gammaRichnessEstimator will produce a set of estimates for the species
% richness of a region given a set of abundance vectors, one for each
% quadrat using the Gamma method as presented in the paper "Community
% Structure to Improve Diversity Estimation" by Skipton Woolley, myself
% (Owen Stewart) and Michael Bode
% this method involves taking the abundance vectors, clustering them into
% communities using kmeans, estimating the species richness of each
% community separately and then accumulating results while correcting for
% overlap between communities
% it will assume the dataset has at least 2 underlying communities present

% note - all inputs apart from abundMat are optional - setting any
% inputs to [] will set them to a default value, as detailed below

% inputs:
% abundMat - a matrix of abundance vectors stored in column format - each
% column is an abundance vector for a different quadrat
% nRuns - the number of times to repeat the Gamma method - i.e. number of
% times to determine the optimal number of clusters, cluster with this
% number of clusters and then estimate species richness (default = 40)
% logInd - if equal to "log", abundMat will be log-transformed before
% clustering occurs using log(abundMat + 1) - this can reduce the impact of
% common species and increase the impact of uncommon species on the
% clustering and is recommended (default = "log")
% kVec - specificies the range of k values (number of clusters) to trial in
% vector form (default = [1 ... max(number of quadrats / 3, 25)]
% replicates - the number of times to re - run the k means clustering on a
% given iteration before choosing the best outcome (default = 20)
% metric - specifies the metric used to determine the optimal k value,
% where "Dunn's average" is the default recommended method used in the
% paper
% "Silhouette" will use the average silhouette score
% "Dunn's average" will use Dunn's index with average pairwise distances
% "Dunn's minmax" will use Dunn's index with maximum and minimum pairwise
% distances
% "AIC" or "BIC" will use the AIC or BIC assuming each quadrat's raw
% abundance data (without being log transformed) are the results given by a
% draw from one of k multinomial distributions, one for each cluster,
% where the mutlinomial probabilities are the average proportions of each
% species over the quadrats assigned to that cluster
% "Dunn's average centroid" will use Dunn's index with average distances
% to the centroid
% "Dunn's average rep" will use the Dunn's average method, but repeat it 13
% times and choose the highest estimate
% "gap", "Davies-Bouldin" and "Calinski-Harabasz" will use the gap
% statistic, Davies-Bouldin index and Calinski-Harabsz index respectively
% (default = "Dunn's average")

% outputs:
% estimateVec - a row vector of outputs from the Gamma estimate method
% optKVec - a row vector of the optimal k values found during each
% iteration

% determine the number of quadrats
nQuads = size(abundMat, 2);

% create an incidence matrix
incMat = abundMat > 0;

% if nRuns has not been specified, set as default to 40
if nargin < 2 || isempty(nRuns)
    nRuns = 40;
end

% set logInd to "log" if it has not been specified
if nargin < 3 || isempty(logInd)
    logInd = "log";
end

% check kVec if it has been specified by user, and if not create it
if nargin < 4 || isempty(kVec)
    
    % want kVec to go up to 25 if the number of quadrats is above 75,
    % otherwise want to go up to 1/3 of the number of quadrats
    if nQuads > 75
        kVec = 2:25;
    else
        kVec = 2:(ceil(nQuads/3));
    end
    
else

    % if kVec has been specified by user, remove any values over 25 for
    % computational complexity reasons
    kVec = kVec(kVec <= 25);
    
end

% if replicates have not been specified, set to 20
if nargin < 5 || isempty(replicates)
    replicates = 20;
end

% set the metric to "Dunn's average" if unspecified or invalid
if nargin < 6 || isempty(metric)
    metric = "Dunn's average";
elseif ~ismember(metric, ["Silhouette", "Dunn's average", "Dunn's minmax", "Dunn's average centroid", "AIC", "BIC", "Dunn's average rep", "gap", "Davies-Bouldin", "Calinski-Harabasz"]) % ---------------------------------UPDATE THIS LATER-----------------------
    warning("Unrecognised metric input. Repeated Dunn's Index using average pairwise distances is now being used.");
    metric = "Dunn's average";
end

% initialise results storage
estimateVec = zeros(1, nRuns);
optKVec = zeros(1, nRuns);

for i = 1:nRuns

    % determine and store optimal k value
    optKVec(i) = determineKKMeans(abundMat, kVec, logInd, metric);
    
    % cluster with the optimal k value, once again checking if a log
    % transform is necessary
    if strcmp(logInd, "log")
        idx = kmeans(log(abundMat' + 1), optKVec(i), 'Replicates', replicates);
    else
        idx = kmeans(abundMat', optKVec(i), 'Replicates', replicates);
    end

    % apply the gamma estimator here
    bCell = createbCell(incMat, idx);
    alphaHatVec = Chao2Clusters(incMat, abundMat, idx);
    estimateVec(i) = gammaEstimator(alphaHatVec, bCell);
    
end

end
