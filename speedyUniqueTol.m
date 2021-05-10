function [u,ia]=speedyUniqueTol(vals,uniquetolval)

% sv=sum(vals,2);
uvals=round(vals/uniquetolval)*uniquetolval;
[u_,ia_]=unique(uvals,'rows');
u=vals(ia_,:);
ia=cell(size(ia_));
%only needed for MEX
% for i=1:size(ia,1)
%    ia{i}=i;
% end
parfor i=1:size(u_,1)
    t=u_(i,:)==uvals;
    ia{i}=find(all(t,'all'));
end