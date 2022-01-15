function V = invCategorical(X,catValues,catIndex)
%restore categorical variables based on index

V = categorical(zeros(size(X)));
for k=1:length(catIndex)
    t = X==catIndex(k);
    V(t) = catValues(k);
end

end