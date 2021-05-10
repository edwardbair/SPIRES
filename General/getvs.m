function vns = getvs(v,number)
%Example vns=getvs(whos,1)%Gets the largest variable
%Example vns=getvs(whos,2)%Gets the 2 largest variables
%INPUT
% use the function 'whos' here
% number - number of variables desired
%OUTPUT
% vns - variable and size sorted by size in GB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Karl Rittger
% UCSB Snow Hydrology
% April 12, 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(v)
   vn{i,1}=v(i).name;
   vs{i,1}=v(i).bytes/ 1073741824;
end

% Sort ascending, flip to descending
[~,indSort] = sortrows(vs);
indSort=flipud(indSort);

for i=1:number
    vns{i,1}=vn{indSort(i)};
    vns{i,2}=[num2str(vs{indSort(i)}) ' GB'];
    vns{i,3}=vs{indSort(i)};
end