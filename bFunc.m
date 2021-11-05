function sum = bFunc(bVec, j, n, bCell, N, testVec)
% bFunc is a recursive function which aims to calculate all beta groupings
% associated with a given community, and then add or subtract them based on
% the equation for the gamma diversity of a region
% it only needs to be called explicitly once for each community - after
% this initial call, it recursively calls itself, moving in a depth first
% search ordering to calculate and sum all relevant beta groupings for a 
% given community

% inputs:
% bVec - the vector representing the current beta intersection - for the
% initial call for community z, bVec should simply be an indexed abundance
% vector for community z
% bVec should be in an indexed format - i.e. rather than a binary abundance
% vector, it should instead be an index of all nonzero elements of its
% corresponding binary incidence vector (i.e. abundVec = [1 0 1 0 0 1]
% needs to be converted to bVec = [1 3 6], which can be done using the find
% command in matlab)
% j - the highest community represented in bVec - for an initial call on a
% community z, we have j = z
% n - the current 'level' - i.e. the number of communities in each group
% which will be calculated explicitly INSIDE this current call - hence bVec
% contains only n - 1 communities, and groups of n communities will be
% calculated in this call through an intersection with bVec 
% groups of > n communities will also be calculated implicitly, using 
% recursive calls whose sums are cumulated in the sum2 variable
% for the initial call, n should be set to 2
% bCell - a cell array which stores the incidence vectors for each
% community in an indexed format identical to that detailed for bVec above
% N - the total number of communities in the region
% testVec - keeps track of the communities represented in bVec - this is
% primarily used for testing, and is an optional parameter

% output:
% sum - the contribution to the overall beta sum provided by the current
% bVec, and all its children

% NOTE: this code does NOT scale the betas by the alpha richness of the
% lowest community represented in bVec - as such, when called initially on
% a community z, the resulting sum should be later divided by alpha_z when
% used in the equation for gamma richness

if isempty(bVec)
    
    % if the beta vector is empty, simply set the sum to 0 and exit the
    % method without making any recursive calls, as all other future
    % intersects with this beta vector would also be empty
    sum = 0;
    
elseif j == N
    
    % if the largest value in the array is equal to the total number of
    % communities, then there are no other possible joins which can be
    % calculated, so set the sum to 0 and exit the method
    sum = 0;
    
else
    
    % otherwise, combinations can be made
    % create two sum variables:
    % sum1 will be used to add up the beta
    % values at the current level, i.e. those which contain n
    % communities
    % sum2 will be used to add up all the beta values at levels n + 1 and
    % beyond, gathered through recursive calls
    sum1 = 0;
    sum2 = 0;
    
    % loop over the possible communities which could be intersected with
    % the current bVec, i.e. any communities higher than j
    for k = (j + 1):N
        
        % determine the intersection between bVec and community k
        bVec2 = intersect(bVec, bCell{k});
        
        % add the number of species common to both to sum1
        sum1 = length(bVec2) + sum1;
        
        % for testing purposes:
        if nargin == 6
            fprintf([num2str([testVec k]) '    length: ' num2str(length(bVec2)) '\n'])
        end
        
        % recursively call for the next level (n + 1 communities), adding 
        % the returned value to sum2
        if nargin == 6
            sum2 = sum2 + bFunc(bVec2, k, n + 1, bCell, N, [testVec k]);
        else
            sum2 = sum2 + bFunc(bVec2, k, n + 1, bCell, N);
        end
    end
    
    % calculate the total sum for this call, applying the correct sign to
    % the betas calculated at the current level
    sum = ((-1)^(n-1)) * sum1 + sum2;
    
end

end