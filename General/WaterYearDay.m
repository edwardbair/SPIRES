function [WYear, DoWY] = WaterYearDay(DatesNum,varargin)
% [WYear, DoWY] = WaterYearDay(DatesNum [,startMonth])
% convert datenum to WY and DOWY
%
% input
%   DatesNum - vector of MATLAB datenum values, optionally can be datetime
%       values
% optional input
%   startMonth - first month of the water year, default 10 (October)
% output, vectors of equal length
%   WYear - water year
%   DoWY - day of water year (starts 1st day of startMonth)

% parse input
p = inputParser;
defaultMonth = 10; % October
addRequired(p,'DatesNum',@(x) isnumeric(x) || isdatetime(x));
addOptional(p,'startMonth',defaultMonth,@isnumeric);
parse(p,DatesNum,varargin{:});

if isnumeric(p.Results.DatesNum)
    DatesNum = p.Results.DatesNum;
else
    DatesNum = datenum(p.Results.DatesNum);
end

% make input a column vector if it's a row vector
multidim = false;
doTranspose = false;
if numel(DatesNum)>1
    if isrow(DatesNum)
        DatesNum = DatesNum.';
        doTranspose = true;
    else % matrix or cube
        nDates = size(DatesNum);
        DatesNum = reshape(DatesNum,numel(DatesNum),1);
        multidim = true;
    end
end

% Get date vector for all dates(y,m,d,H,M,S)
DatesVec=datevec(DatesNum);
% set water year output vector
WYear=DatesVec(:,1);
t=DatesVec(:,2)>=p.Results.startMonth;
if p.Results.startMonth~=1 % don't add if WY = calendar year
    WYear(t)=WYear(t)+1;
end
% initialize zero date as the last day of the previous water year (e.g. 30
% Sept if water year starts in October)
dvinit=zeros(size(DatesVec));
dvinit(:,1)=WYear-1;
if p.Results.startMonth==1
    dvinit(:,2) = 12;
else
    dvinit(:,2) = p.Results.startMonth-1;
end
dvinit(:,3) = eomday(dvinit(:,1),dvinit(:,2));

% Calculate the Day of Water Year.
% When the input DatesNum value has a fractional part, so will the Day of
% Water Year. To prevent this, we use floor().
DoWY=floor(DatesNum-datenum(dvinit));

% if multi-dimensional or originally a row vector, reshape
if multidim
    WYear = reshape(WYear,nDates);
    DoWY = reshape(DoWY,nDates);
elseif doTranspose
    WYear = WYear.';
    DoWY = DoWY.';
end

end