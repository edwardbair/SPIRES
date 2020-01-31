function mstruct = createUTMmstruct(p1,geoid,varargin)
% mstruct = createUTMmstruct(p1,geoid [,units])
%Function to Create UTM projection struct
%
%INPUT
%   p1 - vector length of 2 with lat,long in desired UTM zone
%   geoid - character string with desired geoid
%           Note: MATLAB help alamanac lists ellipsoids
%OPTIONAL INPUT
%   units - character string with desired units, meters is default
%
%OUTPUT
%   mstruct - Map projection structure from given point and geoid
%

%Calculates UTM zone
z1 = utmzone(p1);

%default utm mstruct
mstruct=defaultm('utm');

%set zone point p1
mstruct.zone=z1;

%set geoid
if nargin<3
    units = 'meters';
else
    units = varargin{1};
end
mstruct.geoid=almanac('earth',geoid,units);

%Set appropriate defaults based on parameters in mstruct
mstruct = defaultm(mstruct);
end