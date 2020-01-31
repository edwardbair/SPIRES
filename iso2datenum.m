function matdate = iso2datenum( isodate )
% matdate = iso2datenum( isodate )
%converts ISO dates in form YYYYMMDD (8 digits) or YYYYDDD (7 digits) to
%MATLAB datenums
%
% input isodate must be numeric, not a character string, but can be a vector
% only works for modern dates (i.e. 4 digits in the year)
% works with a mixture of 7- and 8-digit ISO dates

assert(isnumeric(isodate),...
    'input isodate must be numeric, not character string')
isodate = round(isodate); % just in case decimal date entered

% make a column vector to use vectors of year, month, day
if length(isodate)>1 && isrow(isodate)
    % will switch back to row vector at end
    isodate = isodate';
    switchrows = true;
else
    switchrows = false;
end

matdate = zeros(size(isodate));

% some isodates are 7 digits (YYYYDDD), others 8 (YYYYMMDD)
digits = 1+floor(log10(isodate));
d7 = digits==7;
d8 = digits==8;
badDates = ~(d7 | d8);

if nnz(d7) % form YYYYDDD
    y = floor(isodate(d7)/1000);
    d = mod(isodate(d7),1000);
    matdate(d7) = datenum(y-1,12,31)+d;
end

if nnz(d8) % form YYYYMMDD
    y = floor(isodate(d8)/10000);
    m = floor(mod(isodate(d8),10000)/100);
    d = mod(isodate(d8),100);
    matdate(d8) = datenum(y,m,d);
end

if nnz(badDates)
    warning('some bad dates set to NaN')
    warning('%d ',isodate(badDates));
    matdate(badDates) = NaN;
end

if switchrows
    matdate = matdate';
end
end