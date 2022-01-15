function mstruct = createUTMmstruct(p1,varargin)
% mstruct = createUTMmstruct(p1,geoid [,units])
%Function to Create UTM projection struct
%
%INPUT
%   p1 - vector length of 2 with lat,long in desired UTM zone
%Optional input
%   geoid - default 'wgs84'
%
%OUTPUT
%   mstruct - Map projection structure from given point and geoid
%

p = inputParser;
addRequired(p,'p1',@(x) isnumeric(x) && length(x)==2)
addOptional(p,'geoid',7030,@(x) isnumeric(x) || ischar(x));
parse(p,p1,varargin{:});
geoid = p.Results.geoid;

%Calculates UTM zone
z1 = utmzone(p1);

%default utm mstruct
mstruct=defaultm('utm');

%set zone point p1
mstruct.zone=z1;

E = referenceEllipsoid(geoid);
mstruct.geoid = [E.SemimajorAxis E.Eccentricity];

%Set appropriate defaults based on parameters in mstruct
mstruct = defaultm(mstruct);
end