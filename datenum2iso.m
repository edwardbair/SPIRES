function isodate = datenum2iso( matdate, digits )
% isodate = datenum2iso( matdate, digits )
%converts MATLAT datenum(s) to isodate(s), either YYYYDDD or YYYYMMDD
%
% input
%   matdate - MATLAB datenum(s)
%   digits - 7 for YYYYDDD, 8 for YYYYMMDD

assert(digits==7 || digits==8,'input argument digits must be 7 or 8')

% process as column vector, but switch at end if input is row vector
matdate = floor(matdate);
if length(matdate)>1 && isrow(matdate)
    switchrows = true;
    matdate = matdate';
else
    switchrows = false;
end

% year month day of input
X = datevec(matdate);
assert(~any(X(:,1)<1000),...
    'function is restricted to dates on or after 1000-01-01')
switch digits
    case 7 % form YYYYDDD
        isodate = X(:,1)*1000+(matdate-datenum(X(:,1)-1,12,31));
    case 8 % form YYYYMMDD
        isodate = X(:,1)*10000+X(:,2)*100+X(:,3);
end

if switchrows
    isodate = isodate';
end

end