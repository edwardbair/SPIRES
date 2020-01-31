% snagged from MATLAB mle routine
function expanded = expandInput(input,freq)
%EXPANDDATA Expand out an input vector using element frequencies.
assert(isequal(size(input),size(freq)),...
    'expandInput:InputSizeMismatch',...
        'Input argument sizes must match.');
i = cumsum(freq);
j = zeros(1, i(end));
j(i(1:end-1)+1) = 1;
j(1) = 1;
expanded = input(cumsum(j));
end
