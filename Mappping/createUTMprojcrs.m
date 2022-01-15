function proj = createUTMprojcrs(lat,lon)
% proj = createUTMprojcrs(lat,lon)
%Function to Create UTM projection coordinate reference system
%
%INPUT
%   lat,lon - in degrees, lat + N, lon + E
%
%OUTPUT
%   proj - projcrs object
%

% good for versions R2020b and later
assert(~verLessThan('map','5.0'),...
    'function applies only to MATLAB versions R2020b (Mapping Toolbox 5.0) and later')

%Calculates UTM zone
z1 = utmzone(lat,lon);
zone = str2double(z1(1:end-1));
if lat>=0
    id = 326*100+zone;
else
    id = 327*100+zone;
end

% return the projcrs
proj = projcrs(id);

end