function [ allDays ] = allDaysWaterYear( wy, varargin )
% [ allDays ] = allDaysWaterYear( wy, [startmo, startday] )
%   returns vector of all days in specified water year
%
% Input
%   wy - water year, NNNN, with default start date 10/1/(wy-1)
%
% Optional input
%   2 arguments, month and day to start water year, default 10, 1
%
% Output
%   allDays - vector of days (MATLAB datenum format) in this water year

% process optional arguments
startmo = 10;
startday = 1;
optarg = size(varargin,2);
if optarg
    startmo = varargin{1};
    startday = varargin{2};
    assert(startmo>=1 && startmo<=12,'start month must be between 1 & 12')
    assert(startday>=1,'start day cannot be zero or negative')
end

nMonths = 12;
if startday==1
    endmo = mod(startmo+nMonths-1,nMonths);
    endday = eomday(wy,endmo);
else
    if startday>eomday(wy-1,startmo)
        warning('%d-%02d-%02d not possible start day, reset to %02d/01',...
            wy-1,startmo,startday,startmo);
        startday = 1;
        endmo = mod(startmo+nMonths-1,nMonths);
        endday = eomday(wy,endmo);
    else
        endday = startday-1;
        endmo = startmo;
    end
end

allDays = (datenum([wy-1 startmo startday]) : datenum([wy endmo endday]))';

end