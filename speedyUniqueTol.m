function [u,ia]=speedyUniqueTol(vals,uniquetolval)

sv=sum(vals,2);
uvals=round(sv/uniquetolval)*uniquetolval;
[u_,ia_]=unique(uvals);
u=vals(ia_,:);
ia=cell(size(ia_));
%only needed for MEX
% for i=1:size(ia,1)
%    ia{i}=i;
% end
for i=1:size(u_,1)
    ia{i}=find(u_(i)==uvals);
end