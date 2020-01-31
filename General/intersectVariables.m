function [ vList ] = intersectVariables(candidates,masterList)
%which variables are common to both lists
% [ vList ] = intersectVariables( candidates,masterList )
%
% Input
%   candidates - cell list of candidate variables that might be in the
%   masterList, also a cell
%
% Output
%   vList - list of candidates that are in the masterList

vList = {};
for k=1:length(candidates)
    list = strcmpi(candidates{k},masterList);
    if nnz(list)>0
        vList = cat(2,vList,candidates{k});
    end
end
end