function [T,hdr,S]=readSVC(fname)
% [T,hdr,S]=readSVC(fname)
%read in SVC spectra
%input: fname - filename
%output: T - spectra table
%hdr - file header info
%S - struct w matdate,lat & lon of spectra
N=27;
out=readmatrix(fname,...
    'NumHeaderLines',N,...
    'FileType','text',...
    'ExpectedNumVariables',4);
fid=fopen(fname);
hdr=textscan(fid,'%s');

runits='W m^-2 sr^-1 nm^-1';
T=table(out(:,1),out(:,2),out(:,3),out(:,4)./100,...
    'VariableNames',{'wavelength','Irradiance',...
    'RadianceTarget','ReflectanceTarget'});
T.Properties.VariableUnits=...
    {'nm',runits,runits,'unitless'};
if length(hdr{1})==4099
    offset=0;
else
    offset=1;
end
S.matdateUTC=datenum([hdr{1}{78},' ' hdr{1}{92+offset}],...
    'mm/dd/yyyy HHMMSS.FFF');
lat=hdr{1}{89};
S.lat=str2double(lat(1:2))+str2double(lat(3:4))/60+...
    str2double(lat(6:9))/10000/(60*60);
lon=hdr{1}{86-offset};
if lat(end)=='S'
    S.lat=-S.lat;
end
S.lon=str2double(lon(1:3))+str2double(lon(4:5))/60+...
    str2double(lon(7:10))/10000/(60*60);
if lon(end)=='W'
    S.lon=-S.lon;
end

% add datetime field
S.datetime = datetime(S.matdateUTC,'ConvertFrom','datenum','TimeZone','Etc/GMT');
end
