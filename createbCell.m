function bCell = createbCell(incMat, idx)
% createbCell creates a cell array which holds a vector containing the
% species codes present in each community, for each community

% inputs:
% incMat - a matrix of incidence vectors stored in column format - each
% column is an incidence vector for a different partition (i.e. quadrat)
% each incidence vector is stored in the same order as the partArr - i.e.
% the incidence vector for the 3rd partition will be stored in the 3rd
% column of incMat
% idx - the index vector, which specifies the cluster/community each
% incidence vector belongs to

% output:
% bCell - a cell array which stores incidence vectors for each community in
% an indexed format - i.e. if a community contains species 2, 3, and 5 out
% of 5 known species, bCell will hold the vector [2 3 5] for that
% community, rather than [0 1 1 0 1]
% bCell will be ordered community by community, i.e. bCell{3} will
% provide the indexed incidence vector for community 3
% the indexed incidence vectors will be stored as column vectors, the only
% time bCell is used is with the intersect command in bFunc, for which the
% orientation actually doesn't matter anyway

% determine the number of communities
k = max(idx);

% initialise bCell
bCell = cell(1, k);

% loop over each community
for i = 1:k
    
    % create a mask indicating which incidence vectors belong to community
    % i
    mask = (idx == i);
    
    % sum the incidence vectors corresponding to the community i to create
    % a vector holding repeated incidence data
    repIncVec = sum(incMat(:, mask), 2);
    
    % convert this repeated incidence vector into a single incidence vector
    % for the entire community
    incVec = (repIncVec > 0);
    
    % use the find command to index this incidence vector in the form
    % detailed above, and then store in the corresponding cell
    bCell{i} = find(incVec);
    
end

end

