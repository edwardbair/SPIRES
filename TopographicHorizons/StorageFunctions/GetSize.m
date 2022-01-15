function totSize = GetSize(this)
props = properties(this);
if isempty(props)
    totSize = size(this);
else
    totSize = 0;
    for ii=1:length(props)
        currentProperty = this.(props{ii}); %#ok<NASGU>
        s = whos('currentProperty');
        totSize = totSize + s.bytes;
    end
    totSize = [totSize 1];
end
end