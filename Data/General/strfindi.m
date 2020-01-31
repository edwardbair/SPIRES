function k = strfindi( str, pattern )
% k = strfindi( str, pattern )
%like MATLAB function strfind() but ignores case
%
% calls the MATLAB function strfind() but first converts both the string
% str and the pattern to lower case, thus making this a case-insensitive
% implementation of strfind(), which is explicitly case-sensitive

% strfind() returns k as a cell vector if an input is a cell
str = lower(str);
pattern = char(lower(pattern));

k = strfind(str,pattern);

end