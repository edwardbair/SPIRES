function [T,hdr,S]=readSVCjd(fname,varargin)
% [T,hdr,S]=readSVC(fname [,TimeZone])
%read in SVC spectra
%input: fname - filename
%optional input - TimeZone, default 'Etc/GMT+8'
%output: T - spectra table
%hdr - file header info
%S - struct w matdate,lat & lon of spectra

p = inputParser;
addRequired(p,'fname',@ischar);
addOptional(p,'TimeZone','Etc/GMT+8',@ischar);
parse(p,fname,varargin{:});
TimeZone = p.Results.TimeZone;

N=27;
out=readmatrix(fname,...
    'NumHeaderLines',N,...
    'FileType','text',...
    'ExpectedNumVariables',4);
fid=fopen(fname);
hx=textscan(fid,'%s');
hdr = hx{1};

runits='W m^-2 sr^-1 nm^-1';
T=table(out(:,1),out(:,2),out(:,3),out(:,4)./100,...
    'VariableNames',{'wavelength','Irradiance',...
    'RadianceTarget','ReflectanceTarget'});
T.Properties.VariableUnits=...
    {'nm',runits,runits,''};

% search header for time
ktime = strfind(hdr,'time=');
kx = false(size(ktime));
for n=1:length(kx)
    if ~isempty(ktime{n})
        kx(n) = true;
    end
end
kidx = find(kx);
kidx = kidx(2); % time we want is after 'scan time=' and before 'gpstime='
S.datetime = datetime([hdr{kidx+1} ' ' hdr{kidx+2} ' ' hdr{kidx+3}],...
    'InputFormat','MM/dd/yyyy hh:mm:ss a,','TimeZone',TimeZone);
dhold = S.datetime;
dhold.TimeZone = 'Etc/GMT';
S.matdateUTC = datenum(dhold);

% latitude and longitude
klon = strfind(hdr,'longitude=');
kx = false(size(klon));
for n=1:length(kx)
    if ~isempty(klon{n})
        kx(n) = true;
    end
end
kidx = find(kx);
lon = hdr{kidx+1};
xlon = str2double(lon(1:end-1));
S.lon = floor(xlon/100)+floor(mod(xlon,100))/60+(xlon-floor(xlon))/60;
if lon(end)=='W'
    S.lon = -S.lon;
end
klat = strfind(hdr,'latitude=');
kx = false(size(klat));
for n=1:length(kx)
    if ~isempty(klat{n})
        kx(n) = true;
    end
end
kidx = find(kx);
lat = hdr{kidx+1};
xlat = str2double(lat(1:end-1));
S.lat = floor(xlat/100)+floor(mod(xlat,100))/60+(xlat-floor(xlat))/60;
if lat(end)=='S'
    S.lat = -S.lat;
end

% change missing lat/lon from NaN to empty
if isnan(S.lat)
    S.lat = [];
end
if isnan(S.lon)
    S.lon = [];
end
end